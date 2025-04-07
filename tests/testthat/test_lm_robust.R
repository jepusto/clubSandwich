context("lm_robust objects")
set.seed(20190513)

skip_if_not_installed("estimatr")

library(estimatr)

data(mtcars)

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)

# =============== vcovCR ===============

test_that("vcovCR works", {
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2") # works
  vcov_lmr <- vcovCR(lm_robust, ChickWeight$Chick, "CR2") # doesn't work
  
  expect_equal(vcov_lm, vcov_lmr)
})

# =============== model_matrix() ===============

test_that("model_matrix() works", {
  mm_fit <- model_matrix(lm_fit) # works
  mm_rob <- model_matrix(lm_rob) # doesn't work
  
  expect_equal(mm_fit, mm_rob)
})

# =============== residuals_CS() ===============

test_that("residuals_CS() works", {
  rcs_fit <- residuals_CS(lm_fit) # works
  rcs_rob <- residuals_CS(lm_rob) # gives null
  
  expect_equal(rcs_fit, rcs_rob)
})

# =============== coef() ===============

test_that("coef() works", {
  coef_fit <- coef(lm_fit) # works
  coef_rob <- coef(lm_rob) # works?
  
  expect_equal(coef_fit, coef_rob) # works!!
})

# =============== nobs() ===============

test_that("nobs() works", {
  nobs_fit <- nobs(lm_fit) # works
  nobs_rob <- nobs(lm_rob) # works
  
  expect_equal(nobs_fit, nobs_rob) # works!!!
})


# =============== targetVariance() ===============

test_that("targetVariance() works", {
  tV_fit <- targetVariance(lm_fit, ChickWeight$Chick) # works
  tV_rob <- targetVariance(lm_rob, ChickWeight$Chick) # also works
  
  expect_equal(tV_fit, tV_rob) # works!
})

# =============== weightMatrix() ===============

test_that("weightMatrix() works", {
  wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick) # works!!
  wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick) # works!
  
  expect_equal(wM_fit, wM_rob) # works!!!!
})

# =============== v_scale() ===============

test_that("v_scale() works", {
  vs_fit <- v_scale(lm_fit) # works
  vs_rob <- v_scale(lm_rob) # works
  
  expect_equal(vs_fit, vs_rob) # works!!
})
