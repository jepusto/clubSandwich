library(stringr)
library(lubridate)
library(plyr)
library(dplyr)
library(tidyr)
library(plm)
library(nlme)
library(clubSandwich)
rm(list=ls())

# Read in student data files 
stud <- read.table("paper_ClusterRobustTesting/data/TNSTAR/STAR_Students.tab", sep = "\t", header = TRUE)

#-----------------------------------
# Student data 
#------------------------------------

student_dat <- within(stud, {
  female <- ifelse(gender==2, 1, 0)
  minority <- ifelse(race %in% c(1,3), 0, 1)
  white <- 1 - minority
  gk_small <- ifelse(gkclasstype=="1", 1, 0)
  gk_RA <- ifelse(gkclasstype == "3", 1, 0)
  free <- ifelse(gkfreelunch == "1", 1, 0)
  urbanicity <- factor(gksurban, levels = 1:4, c("Inner","Suburb","Rural","Urban"))
  dob <- as.Date(paste(birthyear, birthmonth, birthday, sep = "-"))
  agek <- new_interval(dob, as.Date("1985-09-01")) / duration(num = 1, units = "years")
})

#------------------------------------------
# Outcome measure creation 
#------------------------------------------

# scale outcomes according to percentiles in regular/RA conditions

filter(student_dat, flaggk==1) %>%
  select(stdntid, gkclasstype, gktreadss, gktmathss, gkwordskillss) %>%
  gather("outcome","scale_score",gktreadss, gktmathss, gkwordskillss) %>%
  mutate(outcome = ifelse(outcome=="gktreadss","read",ifelse(outcome=="gktmathss","math","wordskill"))) %>%
  group_by(outcome) %>%
  do(data.frame(stdntid = .$stdntid, test_pct = 100 * ecdf(.$scale_score[.$gkclasstype != 1])(.$scale_score))) %>%
  unite(SID_outcome, stdntid, outcome, remove = FALSE) ->
  student_outcomes

filter(student_dat, flaggk==1) %>%
  select(stdntid, gkschid, white, free, female, agek, gk_small, gk_RA) %>%
  mutate(gkschid = factor(gkschid)) %>%
  left_join(student_outcomes, by = "stdntid") %>%
  mutate(outcome_int = as.integer(factor(outcome))) %>%
  na.omit() ->
  student_dat_outcomes


#------------------------------------------
# fit regressions
#------------------------------------------

# reading

read <- lm(test_pct ~ 0 + gk_small + gk_RA + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="read"))
read_plm <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="read"), index = c("gkschid","stdntid"))
all.equal(coef(read)[1:6], coef(read_plm))

# math

math <- lm(test_pct ~ 0 + gk_small + gk_RA + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="math"))
math_plm <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="math"), index = c("gkschid","stdntid"))
all.equal(coef(math)[1:6], coef(math_plm))

# word recognition

word <- lm(test_pct ~ 0 + gk_small + gk_RA + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="wordskill"))
word_plm <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="wordskill"), index = c("gkschid","stdntid"))
all.equal(coef(word)[1:6], coef(word_plm))

# SUR

SUR <- lm(test_pct ~ 0 + gkschid + outcome + outcome:gk_small + outcome:gk_RA 
          + white + female + free + agek,
          data = student_dat_outcomes)
SUR_plm <- plm(test_pct ~ outcome + outcome:gk_small + outcome:gk_RA + white + female + free + agek,
                data = student_dat_outcomes, index = c("gkschid","SID_outcome"))
all.equal(coef(SUR)[80:91], coef(SUR_plm))

# SUR_gls <- gls(test_pct ~ 0 + gkschid + outcome + outcome:gk_small + outcome:gk_RA
#                + white + female + free + agek, data = student_dat_outcomes, 
#                correlation = corSymm(form = ~ outcome_int | stdntid))
# summary(SUR_gls)$tTable[-(1:79),]

mods <- list(read_lm = read, read_plm = read_plm,
             math_lm = math, math_plm = math_plm,
             word_lm = word, word_plm = word_plm,
             SUR_lm = SUR, SUR_plm = SUR_plm)

#------------------------------------------
# Wald tests
#------------------------------------------

run_Wald_tests <- function(mod) {
  trt_coefs <- which(grepl("gk_", names(coef(mod))))
  if (class(mod) == "lm") {
    clustering_var <- model.frame(mod)$gkschid
    CR0 <- Wald_test(mod, trt_coefs, vcov = "CR0", cluster = clustering_var, test = "Naive-F")
    CR2 <- Wald_test(mod, trt_coefs, vcov = "CR2", cluster = clustering_var, test = c("Naive-F","HTZ")) 
  } else {
    CR0 <- Wald_test(mod, trt_coefs, vcov = "CR0", test = "Naive-F")
    CR2 <- Wald_test(mod, trt_coefs, vcov = "CR2", test = c("Naive-F","HTZ")) 
  }
  CR0$CR <- 0
  CR0$test <- rownames(CR0)
  CR2$CR <- 2
  CR2$test <- rownames(CR2)
  res <- rbind(CR0, CR2)
  class(res) <- "data.frame"
  row.names(res) <- NULL
  select(res, CR, test, Fstat, df, p_val)
}

system.time(F_tests <- ldply(mods, run_Wald_tests, .progress = "text"))
save(F_tests, file = "paper_ClusterRobustTesting/STAR_test_results.Rdata")

trt_coefs <- which(grepl("gk_", names(coef(SUR_gls))))
gls_CR0 <- vcovCR(SUR_gls, type = "CR0")
gls_CR2 <- vcovCR(SUR_gls, type = "CR2")
Wald_test(SUR_gls, trt_coefs, vcov = gls_CR0, test = "Naive-F")
Wald_test(SUR_gls, trt_coefs, vcov = gls_CR2, test = c("Naive-F","HTZ"))


#------------------------
# Test
#------------------------
obj <- SUR
clustering_var <- model.frame(obj)$gkschid
type <- "CR2"
target <- NULL
inverse_var <- FALSE

system.time(CR0 <- vcovCR(obj, cluster = clustering_var, type = "CR0", inverse_var = inverse_var))
system.time(CR1 <- vcovCR(obj, cluster = clustering_var, type = "CR1", inverse_var = inverse_var))
system.time(CR2 <- vcovCR(obj, cluster = clustering_var, type = "CR2", inverse_var = inverse_var))
system.time(CR3 <- vcovCR(obj, cluster = clustering_var, type = "CR3", inverse_var = inverse_var))
system.time(CR4 <- vcovCR(obj, cluster = clustering_var, type = "CR4", inverse_var = inverse_var))
