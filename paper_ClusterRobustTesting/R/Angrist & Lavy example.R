library(plyr)
library(dplyr)
library(tidyr)
library(foreign)
library(plm)
library(sandwich)
library(lmtest)
library(clubSandwich)
devtools::load_all(".")
rm(list=ls())

vcovCRAP <- function(x) vcovHC(x, method = "arellano", type = "HC0")

#-----------------------------------------
# read in and munge data
#-----------------------------------------

filenames <- paste0("paper_ClusterRobustTesting/data/AngristLavy_AERdata/base",c("99","00","01","02"), ".dta")
panel_dat_list <- lapply(filenames, read.dta)
panel_dat <- bind_rows(panel_dat_list)

panel_dat <- within(panel_dat, {
  yr <- factor(year, levels = c(99,0,1,2), labels = 1999:2002)
  school_type <- factor(ifelse(semrel, "Religious", ifelse(semarab, "Arab","Secular")))
  student_id <- paste0(yr, "-", student_id)
  treated2001 <- (yr=="2001") * treated
  sex <- ifelse(boy==1, "Boy","Girl")
  ah4 <- m_ahim >= 4
})

panel_dat <- cbind(panel_dat, model.matrix(~ 0 + yr, data = panel_dat))

qntl <- function(x, groups = 4) {
  qtl <- quantile(x, (0:groups) / groups)
  cut(x, qtl, labels = 1:groups, include.lowest = TRUE, right = FALSE, ordered_result = TRUE)
}

filter(panel_dat, yr!=1999) %>%
  group_by(yr, sex) %>%
  mutate(qrtl = qntl(lagscore, 4),
         half = qntl(lagscore, 2)) %>%
  ungroup() %>%
  arrange(school_id, student_id) ->
  panel_dat

with(panel_dat, table(yr, sex))
with(panel_dat, table(sex, half, yr))

group_by(panel_dat, school_type, treated, school_id) %>%
  summarize(n = n()) %>%
  group_by(school_type, treated) %>%
  summarize(schools = n(), students = sum(n))
#-----------------------------------------
# ATE, panel model
#-----------------------------------------

dat <- filter(panel_dat, sex=="Girl")

ols <- lm(zakaibag ~ factor(school_id) + yr:school_type + 
            educav + educem + ole5 + ah4 + qrtl + treated2001:half, 
          data = dat)

ols_constraints <<- grepl("treated2001", names(coef(ols)))
sum(ols_constraints)
coef_test(ols, vcov="CR1", cluster=dat$school_id, test = "naive-t")[ols_constraints,]
coef_test(ols, vcov="CR2", cluster=dat$school_id, test = c("naive-t","Satterthwaite"))[ols_constraints,]
Wald_test(ols, constraints = ols_constraints, vcov = "CR1", cluster = dat$school_id, test = "Naive-F")
Wald_test(ols, constraints = ols_constraints, vcov = "CR2", cluster = dat$school_id, test = c("Naive-F","HTZ"))


panel <- plm(zakaibag ~ yr2000:school_type + yr2001:school_type + educav + educem + ole5 + ah4 + qrtl + treated2001:half, 
             data = dat, index = c("school_id","student_id"))
summary(panel)
panel_constraints <<- grepl("treated2001", names(coef(panel)))
sum(panel_constraints)

coeftest(panel, vcov. = vcovCRAP)
coef_test(panel, vcov="CR0", test = "z")

waldtest(panel, ~ . - treated2001:half, vcov = vcovCRAP)
Wald_test(panel, panel_constraints, vcov = "CR0", test = "chi-sq")

coef_test(panel, vcov="CR1", test = "naive-t")[panel_constraints,]
coef_test(panel, vcov="CR1", cluster=dat$school_id, test = "naive-t")[panel_constraints,]
coef_test(panel, vcov="CR2", cluster=dat$school_id, test = c("naive-t","Satterthwaite"))[panel_constraints,]
Wald_test(panel, constraints = panel_constraints, vcov = "CR1", cluster = dat$school_id, test = "Naive-F")
Wald_test(panel, constraints = panel_constraints, vcov = "CR2", cluster = dat$school_id, test = c("Naive-F","HTZ"))


#---------------------------------------------------
# school-type by treatment interaction, panel model
#---------------------------------------------------

# test of school-type by treatment interaction

panel_school <- plm(zakaibag ~ yr2000:school_type + yr2001:school_type + educav + educem + ole5 + ah4 + qrtl 
             + treated2001:half + treated2001:half:school_type, 
             data = dat, index = c("school_id","student_id"))
summary(panel_school)
school_constraints <- grepl("school_type.*treated2001", names(coef(panel_school)))

waldtest(panel_school, ~ . - treated2001:half:school_type, vcov = vcovCRAP)
Wald_test(panel_school, school_constraints, vcov = "CR0", test = "chi-sq")

Wald_test(panel_school, constraints = school_constraints, vcov = "CR1", cluster = dat$school_id, test = "Naive-F")
Wald_test(panel_school, constraints = school_constraints, vcov = "CR2", cluster = dat$school_id, test = c("Naive-F","HTZ"))
