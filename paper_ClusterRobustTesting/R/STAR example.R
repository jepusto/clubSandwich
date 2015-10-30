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

read_lm_A <- lm(test_pct ~ 0 + gk_small + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="read"))
read_plm_A <- plm(test_pct ~ gk_small + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="read"), index = c("gkschid","stdntid"))
all.equal(coef(read_lm_A)[1:5], coef(read_plm_A))

read_plm_B <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                  data = filter(student_dat_outcomes, outcome=="read"), index = c("gkschid","stdntid"))

# math

math_lm_A <- lm(test_pct ~ 0 + gk_small + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="math"))
math_plm_A <- plm(test_pct ~ gk_small + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="math"), index = c("gkschid","stdntid"))
all.equal(coef(math_lm_A)[1:5], coef(math_plm_A))

math_plm_B <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                  data = filter(student_dat_outcomes, outcome=="math"), index = c("gkschid","stdntid"))

# word recognition

word_lm_A <- lm(test_pct ~ 0 + gk_small + white + female + free + agek + gkschid,
           data = filter(student_dat_outcomes, outcome=="wordskill"))
word_plm_A <- plm(test_pct ~ gk_small + white + female + free + agek,
                data = filter(student_dat_outcomes, outcome=="wordskill"), index = c("gkschid","stdntid"))
all.equal(coef(word_lm_A)[1:5], coef(word_plm_A))

word_plm_B <- plm(test_pct ~ gk_small + gk_RA + white + female + free + agek,
                  data = filter(student_dat_outcomes, outcome=="wordskill"), index = c("gkschid","stdntid"))

# SUR

SUR_lm_A <- lm(test_pct ~ 0 + gkschid + outcome + outcome:gk_small + white + female + free + agek,
          data = student_dat_outcomes)
SUR_plm_A <- plm(test_pct ~ outcome + outcome:gk_small + white + female + free + agek,
                 data = student_dat_outcomes, index = c("gkschid","SID_outcome"))
all.equal(coef(SUR_lm_A)[80:88], coef(SUR_plm_A))

SUR_plm_B <- plm(test_pct ~ outcome + outcome:gk_small + outcome:gk_RA + white + female + free + agek,
                 data = student_dat_outcomes, index = c("gkschid","SID_outcome"))

mods <- list(read_plm_A = read_plm_A, read_plm_B = read_plm_B,
             math_plm_A = math_plm_A, math_plm_B = math_plm_B,
             word_plm_A = word_plm_A, word_plm_B = word_plm_B,
             SUR_plm_A = SUR_plm_A, SUR_plm_B = SUR_plm_B)

#------------------------------------------
# Wald tests
#------------------------------------------

run_Wald_tests <- function(mod) {
  trt_coefs <- which(grepl("gk_", names(coef(mod))))
  if (all(class(mod) == "lm")) {
    clustering_var <- model.frame(mod)$gkschid
    CR1 <- Wald_test(mod, trt_coefs, vcov = "CR1", cluster = clustering_var, test = "Naive-F")
    CR2 <- Wald_test(mod, trt_coefs, vcov = "CR2", cluster = clustering_var, test = c("Naive-F","HTZ")) 
  } else {
    CR1 <- Wald_test(mod, trt_coefs, vcov = "CR1", test = "Naive-F")
    CR2 <- Wald_test(mod, trt_coefs, vcov = "CR2", test = c("Naive-F","HTZ")) 
  }
  CR1$CR <- 1
  CR1$test <- rownames(CR1)
  CR2$CR <- 2
  CR2$test <- rownames(CR2)
  res <- rbind(CR1, CR2)
  class(res) <- "data.frame"
  row.names(res) <- NULL
  res$q <- length(trt_coefs)
  select(res, q, CR, test, Fstat, df, p_val)
}

system.time(F_tests <- ldply(mods, run_Wald_tests, .progress = "text"))

#------------------------
# Test
#------------------------
obj <- SUR_plm

if ("plm" %in% class(obj)) {
  index <- attr(model.frame(obj),"index")
  clustering_var <- switch(obj$args$effect,
                           individual = index[[1]],
                           time = index[[2]])
} else {
  clustering_var <- model.frame(obj)$gkschid
}
type <- "CR2"
target <- NULL
inverse_var <- FALSE

system.time(CR0 <- vcovCR(obj, cluster = clustering_var, type = "CR0", inverse_var = inverse_var))
system.time(CR1 <- vcovCR(obj, cluster = clustering_var, type = "CR1", inverse_var = inverse_var))
system.time(CR2 <- vcovCR(obj, cluster = clustering_var, type = "CR2", inverse_var = inverse_var))
system.time(CR3 <- vcovCR(obj, cluster = clustering_var, type = "CR3", inverse_var = inverse_var))
system.time(CR4 <- vcovCR(obj, cluster = clustering_var, type = "CR4", inverse_var = inverse_var))


vcov <- CR2
constraints <- which(grepl("gk_", names(coef(obj))))
test <-  c("Naive-F","HTZ")
system.time(walds <- Wald_test(obj, constraints, vcov, test) )
