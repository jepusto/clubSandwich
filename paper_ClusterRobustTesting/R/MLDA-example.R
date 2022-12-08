library(plm)
library(nlme)
library(clubSandwich)

data(MortalityRates)

# subset for deaths in motor vehicle accidents, 1970-1983
MV_deaths <- subset(MortalityRates, cause=="Motor Vehicle" & 
                      year <= 1983 & !is.na(beertaxa), 
                    select = -cause)


#------------------------
# Random effects
#------------------------

RE_fit <- lme(mrate ~ legal + beertaxa + factor(year), data = MV_deaths, random = ~ 1 | state)

RE_CR1 <- Wald_test(RE_fit, constraints = constrain_zero(2), vcov = "CR1", test = "Naive-F")
RE_CR2 <- Wald_test(RE_fit, constraints = constrain_zero(2), vcov = "CR2", test = "HTZ")

#------------------------
# Fixed effects
#------------------------

FE_fit <- plm(mrate ~ legal + beertaxa, data = MV_deaths, 
                             effect = "twoways", index = c("state","year"))
FE_CR1 <- Wald_test(FE_fit, constraints = constrain_zero(1), vcov = "CR1", 
                    cluster = MV_deaths$state, test = "Naive-F")
FE_CR2 <- Wald_test(FE_fit, constraints = constrain_zero(1), vcov = "CR2", 
                    cluster = MV_deaths$state, test = "HTZ")

#------------------------
# Hausmann tests
#------------------------

MV_deaths <- within(MV_deaths, {
  legal_cent <- legal - tapply(legal, state, mean)[factor(state)]
  beertaxa_cent <- beertaxa - tapply(beertaxa, state, mean)[factor(state)]
})

Hausman_fit <- lme(mrate ~ legal + beertaxa + legal_cent + beertaxa_cent + factor(year), 
                   data = MV_deaths, random = ~ 1 | state)
Haus_CR1 <- Wald_test(Hausman_fit, constraints = constrain_zero(4:5), vcov = "CR1", test = "Naive-F")
Haus_CR2 <- Wald_test(Hausman_fit, constraints = constrain_zero(4:5), vcov = "CR2", test = "HTZ")


RE_tests <- bind_rows("Standard" = RE_CR1, "AHT" = RE_CR2, .id = "Test") %>% as.data.frame()
FE_tests <- bind_rows("Standard" = FE_CR1, "AHT" = FE_CR2, .id = "Test") %>% as.data.frame()
Hausman_tests <- bind_rows("Standard" = Haus_CR1, "AHT" = Haus_CR2, .id = "Test") %>% as.data.frame()


MLDA_results <- 
  bind_rows(
    "Random effects" = RE_tests,
    "Fixed effects" = FE_tests,
    "Hausman test" = Hausman_tests, 
    .id = "Hypothesis"
  ) %>%
  select(Hypothesis, Test, "F" = Fstat, df = df_denom, p = p_val) %>%
  mutate(Hypothesis = ifelse(Test=="AHT", NA, Hypothesis))
