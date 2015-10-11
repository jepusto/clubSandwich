library(plyr)
library(dplyr)
library(tidyr)
library(foreign)
library(plm)
library(clubSandwich)
devtools::load_all(".")
rm(list=ls())

#-----------------------------------------
# read in and munge data
#-----------------------------------------

filenames <- paste0("paper_ClusterRobustTesting/data/AngristLavy_AERdata/base",c("99","00","01","02"), ".dta")
base_dat_list <- lapply(filenames, read.dta)
base_dat <- bind_rows(base_dat_list)

base_dat <- within(base_dat, {
  yr <- factor(year, levels = c(99,0,1,2), labels = 1999:2002)
  school_type <- factor(ifelse(semrel, "Religious", ifelse(semarab, "Arab","Secular")))
  student_id <- paste0(yr, "-", student_id)
  treated2001 <- (yr=="2001") * treated
  sex <- ifelse(boy==1, "Boy","Girl")
  ah4 <- m_ahim >= 4
})

group_by(base_dat, yr) %>%
  mutate(qrtl = ntile(lagscore, 4)) %>%
  ungroup() ->
  base_dat

both <- base_dat
both$sex <- "All"

filter(rbind(base_dat, both), yr != "1999") %>%
  group_by(sex) ->
  dat_sex

#-----------------------------------------
# ATE, panel model
#-----------------------------------------

ATE_est <- function(dat) {
  ATE_lm <- lm(zakaibag ~ factor(school_id) + yr:school_type + 
                 educav + educem + ole5 + ah4 + factor(qrtl) + treated2001, 
               data = dat)
  test <- coef_test(ATE_lm, vcov = "CR0", cluster = dat$school_id, test = c("Satterthwaite","saddlepoint"))
  test[rownames(test)=="treated2001",]
}
do(dat_sex, ATE_est(dat = .))


ATE_plm <- plm(zakaibag ~ yr:school_type + educav + educem + ole5 + ah4 + factor(qrtl) + treated2001, 
              data = base_dat, index = c("school_id","student_id"))
coef_test(ATE_plm, vcov = "CR0", test = "Satterthwaite")
obj <- ATE_plm
index <- attr(model.frame(obj),"index")
cluster <- switch(obj$args$effect,
                  individual = index[[1]],
                  time = index[[2]])
target <- NULL
inverse_var <- TRUE

vcov <- vcovCR(obj, type = "CR0")
test <- "Satterthwaite"
Ex_method <- "model"

# test of school-type by treatment interaction
mod6 <- zakaibag ~ 0 + factor(school_id) + yr + educav + educem + ole5 + ah4 + factor(qrtl) + treated2001 + treated2001:school_type
m6 <- lm(mod6, data = filter(dat_sex, sex == "All"))
summary(m6)
do(dat_sex, Wald(dat = ., mod6, "treated2001:"))
