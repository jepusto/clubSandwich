context("ignoring absorbed fixed effects")
set.seed(20190513)

skip_if_not_installed("plm")

library(plm)

data(MortalityRates)
MV_Mortality <- subset(MortalityRates, cause=="Motor Vehicle" & state %in% 1:8)
table(MV_Mortality$state)
MV_Mortality$state_fac <- factor(MV_Mortality$state)
# MV_Mortality$pop <- with(MV_Mortality, 1 + rbinom(nlevels(state_fac), size = 4, prob = 0.5)[state_fac])
summary(MV_Mortality$pop)
MV_Mortality$pop_scale <- with(MV_Mortality, pop / mean(pop))
summary(MV_Mortality$pop_scale)

# model specification

specification <- mrate ~ 0 + legal + beertaxa + beerpercap + winepercap + factor(state)

#-----------------------
# unweighted
#-----------------------

ols_LSDV <- lm(specification, data = MV_Mortality)
ols_within <- plm(update(specification, . ~ . - 0 - factor(state)), data = MV_Mortality, effect = "individual", index = c("state","year"))

test_that("Unweighted lsdv and within estimators are equivalent", {
  lsdv <- coef_test(ols_LSDV, vcov = "CR2", cluster = MV_Mortality$state, coefs = 1:4, p_values = FALSE)
  wthn <- coef_test(ols_within, vcov = "CR2", p_values = FALSE)
  expect_equal(lsdv, wthn)
})

#-----------------------
# iv-weights
#-----------------------

wls_LSDV <- lm(specification, weights = pop_scale, data = MV_Mortality)

MV_Mortality_full <- model.frame(lm(specification, weights = pop_scale, data = MV_Mortality))
U_mat <- model.matrix(update(specification, . ~ . - factor(state)), data = MV_Mortality_full)
T_mat <- model.matrix(~ factor(state), data = MV_Mortality_full)
w <- MV_Mortality_full$"(weights)"
state <- MV_Mortality_full$"factor(state)"
U_absorb <- residuals(stats:::lm.wfit(x = T_mat, y = U_mat, w = w))[,-31]
Y_absorb <- residuals(stats:::lm.wfit(x = T_mat, y = MV_Mortality_full$mrate, w = w))
wls_within <- lm(Y_absorb ~ 0 + U_absorb, weights = w)

test_that("Inverse-variance weighted lsdv and within estimators are equivalent.", {
  lsdv <- coef_test(wls_LSDV, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE, p_values = FALSE)[1:4,]
  wthn <- coef_test(wls_within, vcov = "CR2", cluster = state, inverse_var = TRUE, p_values = FALSE)[1:4,]
  expect_equal(lsdv$beta, wthn$beta, check.attributes = FALSE, tolerance = 10^-5)
  expect_equal(lsdv$SE, wthn$SE, check.attributes = FALSE, tolerance = 10^-2)
  expect_equal(lsdv$tstat, wthn$tstat, check.attributes = FALSE, tolerance = 10^-2)
  expect_equal(lsdv$df_Satt, wthn$df_Satt, check.attributes = FALSE, tolerance = 10^-2)
})

#-----------------------
# p-weights
#-----------------------

test_that("Probability-weighted lsdv and within estimators are not necessarily equivalent.", {
  lsdv <- coef_test(wls_LSDV, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE, coefs = 1:4, p_values = FALSE)
  wthn <- coef_test(wls_within, vcov = "CR2", cluster = state, inverse_var = FALSE, p_values = FALSE)
  expect_equal(lsdv$beta, wthn$beta, check.attributes = FALSE, tolerance = 10^-5)
  expect_equal(lsdv$SE, wthn$SE, check.attributes = FALSE, tolerance = 10^-2)
  expect_equal(lsdv$tstat, wthn$tstat, check.attributes = FALSE, tolerance = 10^-2)
  expect_equal(lsdv$df_Satt, wthn$df_Satt, check.attributes = FALSE, tolerance = 10^-2)
})



