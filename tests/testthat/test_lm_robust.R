context("lm_robust objects")
set.seed(20190513)

skip_if_not_installed("estimatr")

library(estimatr)

# data(mtcars)

data("ChickWeight", package = "datasets")
ChickWeight$wt <- 1 + rpois(nrow(ChickWeight), 3)

lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
wlm_fit <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, data = ChickWeight)
wlm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, data = ChickWeight)

# add these as test cases ^
# add different types as test cases, eg "CR0", "CR1S" = "stata"
# add other options as test cases

# =============== sandwich::bread ===============

test_that("sandwhich::bread works", {
  
  # unweighted tests
  
  bread_lm <- bread(lm_fit)
  bread_rob <- bread(lm_rob)
  
  expect_equal(bread_lm, bread_rob)
  
  bread_wlm <- bread(wlm_fit)
  bread_wrob <- bread(wlm_rob)
  
  expect_equal(bread_wlm, bread_wrob)
})


# =============== model_matrix() ===============

test_that("model_matrix() works", {
  
  # unweighted tests
  
  mm_fit <- model_matrix(lm_fit) 
  mm_rob <- model_matrix(lm_rob)
  
  expect_equal(mm_fit, mm_rob)
  
  # weighted tests
  
  mm_wlm <- model_matrix(wlm_fit) 
  mm_wrob <- model_matrix(wlm_rob)
  
  expect_equal(mm_wlm, mm_wrob)
  
})

# =============== residuals() ===============

test_that("residuals() works", {
  
  # unweighted tests
  
  res_fit <- residuals(lm_fit)
  res_rob <- residuals(lm_rob)
  
  expect_equal(res_fit, res_rob)
  
  # weighted tests
  
  res_wlm <- residuals(wlm_fit)
  res_wrob <- residuals(wlm_rob)
  
  expect_equal(res_wlm, res_wrob)
  
})

# =============== residuals_CS() ===============

test_that("residuals_CS() works", {
  
  # unweighted tests
  
  rcs_fit <- residuals_CS(lm_fit)
  rcs_rob <- residuals_CS(lm_rob)
  
  expect_equal(rcs_fit, rcs_rob)
  
  # weighted tests
  
  rcs_wlm <- residuals_CS(wlm_fit)
  rcs_wrob <- residuals_CS(wlm_rob)
  
  expect_equal(rcs_wlm, rcs_wrob)
})

# =============== coef() ===============

test_that("coef() works", {
  
  # unweighted tests
  
  coef_fit <- coef(lm_fit)
  coef_rob <- coef(lm_rob)
  
  expect_equal(coef_fit, coef_rob)

  # weighted tests
  
  coef_wlm <- coef(wlm_fit)
  coef_wrob <- coef(wlm_rob)
  
  expect_equal(coef_wlm, coef_wrob)
  
})

# =============== nobs() ===============

test_that("nobs() works", {
  
  # unweighted tests
  
  nobs_fit <- nobs(lm_fit)
  nobs_rob <- nobs(lm_rob)
  
  expect_equal(nobs_fit, nobs_rob)
  
  # weighted tests
  
  nobs_wlm <- nobs(wlm_fit)
  nobs_wrob <- nobs(wlm_rob)
  
  expect_equal(nobs_wlm, nobs_wrob)
  
})


# =============== targetVariance() ===============

test_that("targetVariance() works", {
  
  # unweighted tests
  
  tV_fit <- targetVariance(lm_fit, ChickWeight$Chick)
  tV_rob <- targetVariance(lm_rob, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob)
  
  # weighted tests
  
  tV_wlm <- targetVariance(wlm_fit, ChickWeight$Chick)
  tV_wrob <- targetVariance(wlm_rob, ChickWeight$Chick)
  
  expect_equal(tV_wlm, tV_wrob)
})

# =============== weightMatrix() ===============

test_that("weightMatrix() works", {
  
  # unweighted tests
  
  wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick)
  wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob)
  
  # weighted tests
  
  wM_wlm <- weightMatrix(wlm_fit, ChickWeight$Chick)
  wM_wrob <- weightMatrix(wlm_rob, ChickWeight$Chick)
  
  expect_equal(wM_wlm, wM_wrob)
  
})

# =============== v_scale() ===============

test_that("v_scale() works", {
  
  # unweighted tests
  
  vs_fit <- v_scale(lm_fit)
  vs_rob <- v_scale(lm_rob)
  
  expect_equal(vs_fit, vs_rob)
  
  # weighted tests
  
  vs_wlm <- v_scale(wlm_fit)
  vs_wrob <- v_scale(wlm_rob)
  
  expect_equal(vs_wlm, vs_wrob)
  
})


# =============== vcovCR ===============

test_that("vcovCR works", {
  
  # unweighted tests
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR0")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR0")
  
  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR1")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR1")

  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR1p")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR1p")
  
  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR1S")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR1S")
  
  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR2")
  
  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR3")
  vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, "CR3")
  
  expect_equal(vcov_lm, vcov_lmr)
  expect_equal(lm_rob$vcov, vcov_lmr)
  
  # types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")
  # 
  # for(type in types) {
  #   vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, type)
  #   vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, type)
  #   
  #   expect_equal(vcov_lm, vcov_lmr)
  #   expect_equal(lm_rob$vcov, vcov_lmr)
  # }
  
  # weighted tests

  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR0")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR0")

  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)

  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR1")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR1")
  
  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)
  
  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR1p")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR1p")
  
  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)
  
  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR1S")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR1S")
  
  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)
  
  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR2")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR2")
  
  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)
  
  vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR3")
  vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, "CR3")
  
  expect_equal(vcov_wlm, vcov_wlmr)
  expect_equal(wlm_rob$vcov, vcov_wlmr)
  
    
  # for(type in types) {
  #   vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, type)
  #   vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, type)
  #   
  #   expect_equal(vcov_wlm, vcov_wlmr)
  #   expect_equal(wlm_rob$vcov, vcov_wlmr)
  # }
  
})
