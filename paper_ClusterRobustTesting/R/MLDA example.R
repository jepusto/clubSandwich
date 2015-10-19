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

# replicate Table 5.2

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

filter(deaths, agegr=="18-20 yrs" & year <= 1983 & dtype == "MVA" & !is.na(beertaxa)) %>%
  select(-agegr, -dtype) %>%
  group_by(state) %>%
  mutate(legal_s = legal - mean(legal),
         beertaxa_s = beertaxa - mean(beertaxa)) ->
  death_dat

RE1 <- lme(mrate ~ legal + legal_s + factor(year), 
           data = death_dat,
           random = ~ 1 | state)
Wald_test(RE1, constraints = 3, vcov = "CR1", test = "Naive-F")
Wald_test(RE1, constraints = 3, vcov = "CR2", test = c("Naive-F", "HTZ"))

obj <- RE1
constraints <- 3
vcov <- vcovCR(obj, type = "CR1")
test <- "Naive-F"

# no trends, with year FE, controlling for beer taxes
FE_fit_lm <- lm(mrate ~ 0 + legal + beertaxa + factor(state) + factor(year), data = dat)
coef_test(FE_fit, vcov = "CR0", cluster = dat$state, test = "z")[c("legal","beertaxa"),]
FE_fit_plm <- plm(mrate ~ legal + beertaxa, data = dat, 
                  index = c("state","year"),
                  effect = "twoways")
coef_test(FE_fit_plm, vcov = "CR0", cluster = "individual", test = "z")

RE_fit <- lme(mrate ~ legal + beertaxa + legal_s + beertaxa_s + factor(year), 
              data = dat,
              random = ~ 1 | state)
summary(RE_fit)
coef_test(RE_fit, vcov = "CR1", test = "naive-t")[c("legal_s","beertaxa_s"),]
coef_test(RE_fit, vcov = "CR2", test = c("naive-t","Satterthwaite"))[c("legal_s","beertaxa_s"),]
Wald_test(RE_fit, constraints = 4:5, vcov = "CR1", test = "Naive-F")
Wald_test(RE_fit, constraints = 4:5, vcov = "CR2", test = c("Naive-F","HTZ"))
