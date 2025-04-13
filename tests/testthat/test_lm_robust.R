context("lm_robust objects")
set.seed(20190513)

skip_if_not_installed("estimatr")

library(estimatr)

# data(mtcars)

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)

# =============== vcovCR ===============

test_that("vcovCR works", {
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR2") # workS
  
  # should these even be equal? They are not.
  expect_equal(vcov_lm, vcov_lmr)
  # check that vcov_lm and vcov_lmr have an identical structure, but not necessarily equal
  expect_identical(dim(vcov_lm), dim(vcov_lmr)) 
})


# =============== model_matrix() ===============

test_that("model_matrix() works", {
  mm_fit <- model_matrix(lm_fit) 
  mm_rob <- model_matrix(lm_rob) # works as of recently
  
  expect_equal(mm_fit, mm_rob) # works as of recently
})

# =============== residuals() ===============

test_that("residuals() works", {
  res_fit <- residuals(lm_fit)
  res_rob <- residuals(lm_rob) # tweak so that it doesn't need data =
  
  expect_equal(res_fit, res_rob) # doesn't work
})

# =============== residuals_CS() ===============

test_that("residuals_CS() works", {
  rcs_fit <- residuals_CS(lm_fit)
  rcs_rob <- residuals_CS(lm_rob)
  
  expect_equal(rcs_fit, rcs_rob) # doesn't work
})

# =============== coef() ===============

test_that("coef() works", {
  coef_fit <- coef(lm_fit)
  coef_rob <- coef(lm_rob)
  
  expect_equal(coef_fit, coef_rob)
})

# =============== nobs() ===============

test_that("nobs() works", {
  nobs_fit <- nobs(lm_fit)
  nobs_rob <- nobs(lm_rob)
  
  expect_equal(nobs_fit, nobs_rob)
})


# =============== targetVariance() ===============

test_that("targetVariance() works", {
  tV_fit <- targetVariance(lm_fit, ChickWeight$Chick)
  tV_rob <- targetVariance(lm_rob, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob)
})

# =============== weightMatrix() ===============

test_that("weightMatrix() works", {
  wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick)
  wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob)
})

# =============== v_scale() ===============

test_that("v_scale() works", {
  vs_fit <- v_scale(lm_fit)
  vs_rob <- v_scale(lm_rob)
  
  expect_equal(vs_fit, vs_rob)
})

# =============== bread ===============

test_that("bread works.", {

  B_lm <- bread(lm_fit)
  B_lmr <- bread(lm_rob) # workS
  
  # should these even be equal? They are not.
  expect_equal(B_lm, B_lmr)
})
