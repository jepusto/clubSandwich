setwd("paper_ClusterRobustTesting")
library(tidyr)
library(dplyr)
library(foreign)
library(plm)
devtools::load_all()

deaths <- read.dta("data/deaths.dta", convert.factors=FALSE)

filter(deaths, agegr==2 & year <= 1983 & dtype %in% c(1, 2, 3, 6)) %>%
  mutate(death = factor(dtype, labels = c("All deaths","Motor vehicle accidents",
                                          "Suicide","All internal causes"))) %>%
  select(-agegr) %>%
  group_by(death) ->
  death_dat

# replicate Tables 5.2 and 5.3

est <- function(dat, controls = NULL) {
  frml_string <- "mrate ~ 0 + legal"
  if (length(controls) > 0) frml_string <- paste(frml_string, controls, sep = " + ")
  
  notrend_noweight <- lm(as.formula(paste0(frml_string, " + factor(state) + factor(year)")), data = dat)
  A <- coef_test(notrend_noweight, vcov = "CR1s", cluster = dat$state, test = "z")["legal",]
  trend_noweight <- lm(as.formula(paste0(frml_string, " + factor(state) + factor(year) + year:factor(state)")), data = dat)
  B <- coef_test(trend_noweight, vcov = "CR1s", cluster = dat$state, test = "z")["legal",]
  notrend_weight <- lm(as.formula(paste0(frml_string, " + factor(state) + factor(year)")), weights = pop, data = dat)
  C <- coef_test(notrend_weight, vcov = "CR1s", cluster = dat$state, test = "z")["legal",]
  trend_weight <- lm(as.formula(paste0(frml_string, " + factor(state) + factor(year) + year:factor(state)")), weights = pop, data = dat)
  D <- coef_test(trend_weight, vcov = "CR1s", cluster = dat$state, test = "z")["legal",]
  res <- cbind(trend = c("No","Yes","No","Yes"),
               weight = c("No","No","Yes","Yes"),
               rbind(A, B, C, D))
  rownames(res) <- NULL
  res
}

est_5.2 <- do(death_dat, est(.))
est_5.3 <- do(death_dat, est(., controls = "beertaxa"))

#--------------------------------
# Hausman-type tests
#--------------------------------
library(nlme)
deaths <- read.dta("data/deaths.dta")

filter(deaths, agegr=="18-20 yrs" & year <= 1983 & dtype == "MVA" & !is.na(beertaxa)) %>%
  select(-agegr, -dtype) %>%
  group_by(state) %>%
  mutate(legal_s = legal - mean(legal),
         beertaxa_s = beertaxa - mean(beertaxa)) ->
  death_dat

# just policy variable

RE1 <- lme(mrate ~ legal + legal_s + factor(year), 
           data = death_dat,
           random = ~ 1 | state)
summary(RE1)
coef_test(RE1, vcov = "CR1", test = "naive-t")["legal_s",]
coef_test(RE1, vcov = "CR2", test = c("naive-t", "Satterthwaite"))["legal_s",]
Wald_test(RE1, constraints = 3, vcov = "CR1", test = "Naive-F")
Wald_test(RE1, constraints = 3, vcov = "CR2", test = c("Naive-F", "HTZ"))


# controlling for beer taxes

RE2 <- lme(mrate ~ legal + beertaxa + legal_s + beertaxa_s + factor(year), 
              data = death_dat,
              random = ~ 1 | state)
summary(RE2)
coef_test(RE2, vcov = "CR1", test = "naive-t")[c("legal_s","beertaxa_s"),]
coef_test(RE2, vcov = "CR2", test = c("naive-t","Satterthwaite"))[c("legal_s","beertaxa_s"),]
Wald_test(RE2, constraints = 4:5, vcov = "CR1", test = "Naive-F")
Wald_test(RE2, constraints = 4:5, vcov = "CR2", test = c("Naive-F","HTZ"))

#--------------------------------
# Random effects estimates
#--------------------------------

# just policy variable

RE1_fit <- lme(mrate ~ legal + factor(year), 
           data = death_dat,
           random = ~ 1 | state)
summary(RE1_fit)
coef_test(RE1_fit, vcov = "CR1", test = "naive-t")["legal",]
coef_test(RE1_fit, vcov = "CR2", test = c("naive-t", "Satterthwaite"))["legal",]
Wald_test(RE1_fit, constraints = 2, vcov = "CR1", test = "Naive-F")
Wald_test(RE1_fit, constraints = 2, vcov = "CR2", test = c("Naive-F", "HTZ"))

# controlling for beer taxes

RE2_fit <- lme(mrate ~ legal + beertaxa + factor(year), 
           data = death_dat,
           random = ~ 1 | state)
summary(RE2_fit)
coef_test(RE2_fit, vcov = "CR1", test = "naive-t")[c("legal","beertaxa"),]
coef_test(RE2_fit, vcov = "CR2", test = c("naive-t","Satterthwaite"))[c("legal","beertaxa"),]
Wald_test(RE2_fit, constraints = 2, vcov = "CR1", test = "Naive-F")
Wald_test(RE2_fit, constraints = 2, vcov = "CR2", test = c("Naive-F","HTZ"))

#--------------------------------
# fixed effects estimates
#--------------------------------

# just policy variable

FE1_fit <- lm(mrate ~ 0 + legal + factor(year) + factor(state), data = death_dat)
summary(FE1_fit)
coef_test(FE1_fit, vcov = "CR1", cluster = death_dat$state, test = "naive-t")["legal",]
coef_test(FE1_fit, vcov = "CR2", cluster = death_dat$state, test = c("naive-t", "Satterthwaite"))["legal",]
Wald_test(FE1_fit, constraints = 1, vcov = "CR1", cluster = death_dat$state, test = "Naive-F")
Wald_test(FE1_fit, constraints = 1, vcov = "CR2", cluster = death_dat$state, test = c("Naive-F", "HTZ"))

# controlling for beer taxes

FE2_fit <- lm(mrate ~ 0 + legal + beertaxa + factor(year) + factor(state), data = death_dat)
summary(FE2_fit)
coef_test(FE2_fit, vcov = "CR1", cluster = death_dat$state, test = "naive-t")[c("legal","beertaxa"),]
coef_test(FE2_fit, vcov = "CR2", cluster = death_dat$state, test = c("naive-t","Satterthwaite"))[c("legal","beertaxa"),]
Wald_test(FE2_fit, constraints = 1, vcov = "CR1", cluster = death_dat$state, test = "Naive-F")
Wald_test(FE2_fit, constraints = 1, vcov = "CR2", cluster = death_dat$state, test = c("Naive-F","HTZ"))
