context("lm_robust objects")
set.seed(20190513)

skip_if_not_installed("estimatr")

library(estimatr)

# data(mtcars)

data("ChickWeight", package = "datasets")
ChickWeight$wt <- 1 + rpois(nrow(ChickWeight), 3)
ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)

lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                    clusters = Chick)
                    # clusters = factor(ChickWeight$Chick, ordered = FALSE))

wlm_fit <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, data = ChickWeight)
wlm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                     data = ChickWeight, 
                     clusters = Chick)
                     # clusters = factor(ChickWeight$Chick, ordered = FALSE))

# add other lm_robust options as test cases

# =============== sandwich::bread ===============

test_that("sandwhich::bread works", {
  
  # unweighted tests
  
  bread_lm <- bread(lm_fit)
  bread_rob <- bread(lm_rob)
  
  expect_equal(bread_lm, bread_rob)
  
  # weighted tests
  
  bread_wlm <- bread(wlm_fit)
  bread_wrob <- bread(wlm_rob)
  
  expect_equal(bread_wlm, bread_wrob)
  
})

# =============== model.frame() ===============

test_that("model.frame() works", {
  
  # NOTE: These tests are identical from those for model_matrix(), except, they
  # use model.frame()
  
  # unweighted tests
  
  mm_fit <- model.frame(lm_fit) 
  mm_rob <- model.frame(lm_rob)
  
  expect_equal(mm_fit, mm_rob)
  
  # weighted tests
  
  mm_wlm <- model.frame(wlm_fit)
  mm_wrob <- model.frame(wlm_rob)
  
  expect_equal(mm_wlm, mm_wrob)
  
})

# =============== model.matrix() ===============

test_that("model.matrix() works", {
  
  # NOTE: These tests are identical from those for model_matrix(), except, they
  # use model.matrix()
  
  # unweighted tests
  
  mm_fit <- model.matrix(lm_fit) 
  mm_rob <- model.matrix(lm_rob)
  
  expect_equal(mm_fit, mm_rob)
  
  # weighted tests
  
  mm_wlm <- model.matrix(wlm_fit)
  mm_wrob <- model.matrix(wlm_rob)
  
  expect_equal(mm_wlm, mm_wrob)
  
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

# add helper function + apply or revert to for loop

test_that("vcovCR works", {
  
  types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")
  
  # unweighted tests
  
  for(type in types) {
    vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, type)
    vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, type)

    expect_equal(vcov_lm, vcov_lmr)
    if(type == "CR2") {
      class(lm_rob$vcov) <- c("vcovCR", "clubSandwich") # coerce to same class
      expect_equal(lm_rob$vcov, vcov_lm)
      expect_equal(lm_rob$vcov, vcov_lmr)
    }
  }
  
  # weighted tests

  for(type in types) {
    vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, type)
    vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, type)

    expect_equal(vcov_wlm, vcov_wlmr)
    if(type == "CR2") {
      class(wlm_rob$vcov) <- c("vcovCR", "clubSandwich") # coerce to same class
      expect_equal(wlm_rob$vcov, vcov_wlm)
      expect_equal(wlm_rob$vcov, vcov_wlmr)
    }
  }
  
})

test_that("vovCR properly pulls cluster specified for lm_robust", {
  
  # reset ChickWeight
  data("ChickWeight")
  
  # unweighted tests
  
  uw_clust <- vcovCR(lm_rob, ChickWeight$Chick, "CR2")
  uw_no_clust <- vcovCR(lm_rob, type = "CR2")
  uw_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2")
  
  expect_equal(uw_clust, uw_no_clust)
  expect_equal(uw_no_clust, uw_lm)
  
  # create an lm_robust that draws in data differently
  lm_rob_fact <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                           clusters = factor(ChickWeight$Chick, ordered = FALSE))
  # perform vcovCR
  uw_fact_cr <- vcovCR(lm_rob_fact, type = "CR2")
  
  expect_equal(uw_clust, uw_fact_cr)
  
  # weighted tests
  
  w_clust <- vcovCR(wlm_rob, ChickWeight$Chick, "CR2")
  w_no_clust <- vcovCR(wlm_rob, type = "CR2")
  w_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2")
  
  expect_equal(w_clust, w_no_clust)
  expect_equal(w_no_clust, w_lm)
  
  # create an lm_robust that draws in data differently
  lm_rob_fact_w <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                             data = ChickWeight, 
                             clusters = factor(ChickWeight$Chick, ordered = FALSE))
  # perform vcovCR
  w_fact_cr <- vcovCR(lm_rob_fact_w, type = "CR2")
  
  expect_equal(uw_clust, w_fact_cr)
  
})
