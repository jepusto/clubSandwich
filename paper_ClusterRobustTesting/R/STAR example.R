library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
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
  do(data.frame(stdntid = .$stdntid, test_pct = 100 * ecdf(.$scale_score[.$gkclasstype != 1])(.$scale_score))) ->
  student_outcomes

#------------------------------------------
# Seemingly unrelated regression
#------------------------------------------

filter(student_dat, flaggk==1) %>%
  select(stdntid, gkschid, white, free, female, agek, gk_small, gk_RA) %>%
  left_join(student_outcomes, by = "stdntid") %>%
  na.omit() ->
  student_dat_outcomes

SUR <- lm(test_pct ~ 0 + outcome + outcome:gk_small + outcome:gk_RA 
          + white + female + free + agek + as.factor(gkschid),
          data = student_dat_outcomes)
coef(SUR)
vcov_SUR <- vcovCR(SUR, cluster = student_dat_outcomes$gkschid, type = "CR2", inverse_var = TRUE)
Wald_test(SUR, 86:91 , vcov = vcov_SUR, test = "All") 
