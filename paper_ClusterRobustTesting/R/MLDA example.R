library(foreign)
library(nlme)
# setwd("paper_ClusterRobustTesting")
deaths <- read.dta("data/deaths.dta") 

filter(deaths, agegr=="18-20 yrs" & year <= 1983 & dtype == "MVA" & !is.na(beertaxa)) %>%
  select(-agegr, -dtype) %>%
  group_by(state) %>%
  mutate(legal_s = legal - mean(legal),
         beertaxa_s = beertaxa - mean(beertaxa)) ->
  death_dat

# Random effects

RE1_fit <- lme(mrate ~ legal + factor(year), data = death_dat, random = ~ 1 | state)
RE1_CR1 <- Wald_test(RE1_fit, constraints = 2, vcov = "CR1", test = "Naive-F")
RE1_CR2 <- Wald_test(RE1_fit, constraints = 2, vcov = "CR2", test = c("Naive-F", "HTZ"))

RE2_fit <- lme(mrate ~ legal + beertaxa + factor(year), data = death_dat, random = ~ 1 | state)
RE2_CR1 <- Wald_test(RE2_fit, constraints = 2, vcov = "CR1", test = "Naive-F")
RE2_CR2 <- Wald_test(RE2_fit, constraints = 2, vcov = "CR2", test = c("Naive-F","HTZ"))

# Fixed effects

FE1_fit <- lm(mrate ~ 0 + legal + factor(year) + factor(state), data = death_dat)
FE1_CR1 <- Wald_test(FE1_fit, constraints = 1, vcov = "CR1", cluster = death_dat$state, test = "Naive-F")
FE1_CR2 <- Wald_test(FE1_fit, constraints = 1, vcov = "CR2", cluster = death_dat$state, test = c("Naive-F", "HTZ"))

FE2_fit <- lm(mrate ~ 0 + legal + beertaxa + factor(year) + factor(state), data = death_dat)
FE2_CR1 <- Wald_test(FE2_fit, constraints = 1, vcov = "CR1", cluster = death_dat$state, test = "Naive-F")
FE2_CR2 <- Wald_test(FE2_fit, constraints = 1, vcov = "CR2", cluster = death_dat$state, test = c("Naive-F","HTZ"))

X <- model.matrix(FE2_fit)
R <- residuals(lm.fit(X[,-(1:2)], X[,1:2]))
FE2_absorb <- lm(death_dat$mrate ~ 0 + R)
FE2_CR2A <- Wald_test(FE2_absorb, constraints = 1, vcov = "CR2", cluster = death_dat$state, test = c("Naive-F","HTZ"))

# Hausmann tests

RE1_Hausman <- lme(mrate ~ legal + legal_s + factor(year), data = death_dat, random = ~ 1 | state)
Haus1_CR1 <- Wald_test(RE1_Hausman, constraints = 3, vcov = "CR1", test = "Naive-F")
Haus1_CR2 <- Wald_test(RE1_Hausman, constraints = 3, vcov = "CR2", test = c("Naive-F", "HTZ"))

RE2_Hausman <- lme(mrate ~ legal + beertaxa + legal_s + beertaxa_s + factor(year), 
                   data = death_dat, random = ~ 1 | state)
Haus2_CR1 <- Wald_test(RE2_Hausman, constraints = 4:5, vcov = "CR1", test = "Naive-F")
Haus2_CR2 <- Wald_test(RE2_Hausman, constraints = 4:5, vcov = "CR2", test = c("Naive-F","HTZ"))

