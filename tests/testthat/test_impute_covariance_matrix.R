context("impute_covariance_matrix")
set.seed(20190513)

K <- 10
N <- sum(1:K)
dat <- data.frame(study = rep(LETTERS[1:K], 1:K), 
                  yi = rnorm(N), 
                  vi = rchisq(N, df = 2),
                  ti = sample(1:(10 * N), N))

test_that("impute_covariance_matrix error messages and missing argument handling are correct.", {
  expect_error(impute_covariance_matrix(vi = dat$vi, cluster = dat$study))  
  expect_error(impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ti = dat$ti)) 
  expect_error(impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ar1 = 0.8)) 
  
  V1 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6)
  V2 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti)
  V3 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti, ar1 = 0)
  V4 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ti = dat$ti, ar1 = 0.5)
  V5 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.0, ti = dat$ti, ar1 = 0.5)
  
  expect_equal(V1, V2) 
  expect_equal(V1, V3)   
  expect_equal(V4, V5)
})

test_that("impute_covariance_matrix returns correct correlations.", {
  r <- 0.7
  V_single_r <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = r)
  r_list <- rbeta(K, 2, 2)
  V_multiple_r <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = r_list)
  
  check_correlation <- function(M, r) if (nrow(M) > 1) max(abs(cov2cor(M)[lower.tri(M)] - r)) else 0
  check_singles <- sapply(V_single_r, check_correlation, r = r)
  expect_true(all(check_singles < 10^-14))
  check_multiples <- Map(check_correlation, M = V_multiple_r, r = r_list)
  expect_true(all(check_multiples < 10^-14))
  
  dat_scramble <- dat[sample(nrow(dat)),]
  V_mat <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = r)
  expect_equal(dat_scramble$vi, diag(V_mat))
  
  V_resorted <- V_mat[order(dat_scramble$study), order(dat_scramble$study)]
  dat_unscramble <- dat_scramble[order(dat_scramble$study),]
  V_unscramble <- impute_covariance_matrix(vi = dat_unscramble$vi, cluster = dat_unscramble$study, r = r)
  expect_equal(V_resorted, metafor::bldiag(V_unscramble))
})

test_that("impute_covariance_matrix works with unobserved factors.", {
  K <- 10
  N <- sum(1:K)
  dat <- data.frame(study = rep(LETTERS[1:K], 1:K), 
                    yi = rnorm(N), 
                    vi = rchisq(N, df = 2))
  levels(dat$study) <- LETTERS[1:(K + 3)]
  
  r <- 0.7
  V_single_r <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = r)
  r_list <- rbeta(K, 2, 2)
  V_multiple_r <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = r_list)
  
  check_correlation <- function(M, r) if (nrow(M) > 1) max(abs(cov2cor(M)[lower.tri(M)] - r)) else 0
  check_singles <- sapply(V_single_r, check_correlation, r = r)
  expect_true(all(check_singles < 10^-14))
  
  dat_scramble <- dat[sample(nrow(dat)),]
  V_mat <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = r)
  expect_equal(dat_scramble$vi, diag(V_mat))
  
  V_resorted <- V_mat[order(dat_scramble$study), order(dat_scramble$study)]
  dat_unscramble <- dat_scramble[order(dat_scramble$study),]
  V_unscramble <- impute_covariance_matrix(vi = dat_unscramble$vi, cluster = dat_unscramble$study, r = r)
  expect_equal(V_resorted, metafor::bldiag(V_unscramble))
})


test_that("impute_covariance_matrix works with AR1 argument.", {
  
})

test_that("impute_covariance_matrix works with subgroup argument.", {
  
})