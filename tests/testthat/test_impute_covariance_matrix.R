context("impute_covariance_matrix")
set.seed(20190513)

K <- 20
N <- sum(1:K)
dat <- data.frame(
  study = rep(LETTERS[1:K], 1:K), 
  yi = rnorm(N), 
  vi = rchisq(N, df = 2),
  ti = sample(1:(10 * N), N),
  si = sample(letters[1:3], N, replace = TRUE)
)

dat$v_study <- unsplit(tapply(dat$vi, dat$study, mean), dat$study)

dat_scramble <- dat[sample.int(nrow(dat)),]


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
  expect_equal(V_resorted, unblock(V_unscramble))
  
})



test_that("impute_covariance_matrix works with AR1 argument.", {
  
  ar1 <- 0.5
  V_mat <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ti = 1:nrow(dat), ar1 = ar1)
  r_mat <- lapply(V_mat, cov2cor)
  r_sums <- sapply(r_mat, sum)
  check <- sapply(2:K, function(k) k + 2 * sum((1:(k-1)) * ar1^((k-1):1)))
  expect_equal(r_sums[-1], check, check.attributes = FALSE)
  
  V_mat_big <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, 
                                         ti = 1:nrow(dat), ar1 = ar1, return_list = FALSE)
  r_mat_big <- cov2cor(V_mat_big)
  expect_equal(sum(r_mat_big), sum(check) + 1)
  v_sums <- sapply(V_mat, sum)
  expect_equal(sum(v_sums), sum(V_mat_big))
  
})

test_that("impute_covariance_matrix works with subgroup argument.", {
  
  X <- model.matrix(~ si + 0, data = dat)
  X_list <- by(X, dat$study, as.matrix)
  
  check_diag <- function(v, x) {
    XVX <- t(x) %*% v %*% x
    all.equal(sum(XVX), sum(diag(XVX)))
  }
  
  check_all_diag <- function(v_list,x_list) {
    all(mapply(check_diag , v = v_list, x = x_list))
  }
  
  V1 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, subgroup = dat$si)
  V2 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti, subgroup = dat$si)
  V3 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti, ar1 = 0, subgroup = dat$si)
  V4 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ti = dat$ti, ar1 = 0.5, subgroup = dat$si)
  V5 <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.0, ti = dat$ti, ar1 = 0.5, subgroup = dat$si)

  expect_true(check_all_diag(V1, X_list))
  expect_true(check_all_diag(V2, X_list))
  expect_true(check_all_diag(V3, X_list))
  expect_true(check_all_diag(V4, X_list))
  expect_true(check_all_diag(V5, X_list))
  
  X_scramble <- model.matrix(~ si + 0, data = dat_scramble)
  
  V1 <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, subgroup = dat_scramble$si)
  V2 <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti, subgroup = dat_scramble$si)
  V3 <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti, ar1 = 0, subgroup = dat_scramble$si)
  V4 <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, ti = dat_scramble$ti, ar1 = 0.5, subgroup = dat_scramble$si)
  V5 <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.0, ti = dat_scramble$ti, ar1 = 0.5, subgroup = dat_scramble$si)

  expect_true(check_diag(V1, X_scramble))  
  expect_true(check_diag(V2, X_scramble))  
  expect_true(check_diag(V3, X_scramble))  
  expect_true(check_diag(V4, X_scramble))  
  expect_true(check_diag(V5, X_scramble))  
  
})

test_that("impute_covariance_matrix works with smooth argument.", {
  
  V1 <- impute_covariance_matrix(vi = dat$v_study, cluster = dat$study, r = 0.6)
  V2 <- impute_covariance_matrix(vi = dat$v_study, cluster = dat$study, r = 0.6, ti = dat$ti)
  V3 <- impute_covariance_matrix(vi = dat$v_study, cluster = dat$study, r = 0.6, ti = dat$ti, ar1 = 0)
  V4 <- impute_covariance_matrix(vi = dat$v_study, cluster = dat$study, ti = dat$ti, ar1 = 0.5)
  V5 <- impute_covariance_matrix(vi = dat$v_study, cluster = dat$study, r = 0.0, ti = dat$ti, ar1 = 0.5)
  
  V1s <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, smooth_vi = TRUE)
  V2s <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti, smooth_vi = TRUE)
  V3s <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.6, ti = dat$ti, ar1 = 0, smooth_vi = TRUE)
  V4s <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, ti = dat$ti, ar1 = 0.5, smooth_vi = TRUE)
  V5s <- impute_covariance_matrix(vi = dat$vi, cluster = dat$study, r = 0.0, ti = dat$ti, ar1 = 0.5, smooth_vi = TRUE)
  
  expect_equal(V1, V1s)
  expect_equal(V2, V2s)
  expect_equal(V3, V3s)
  expect_equal(V4, V4s)
  expect_equal(V5, V5s)
  
  V1 <- impute_covariance_matrix(vi = dat_scramble$v_study, cluster = dat_scramble$study, r = 0.6)
  V2 <- impute_covariance_matrix(vi = dat_scramble$v_study, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti)
  V3 <- impute_covariance_matrix(vi = dat_scramble$v_study, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti, ar1 = 0)
  V4 <- impute_covariance_matrix(vi = dat_scramble$v_study, cluster = dat_scramble$study, ti = dat_scramble$ti, ar1 = 0.5)
  V5 <- impute_covariance_matrix(vi = dat_scramble$v_study, cluster = dat_scramble$study, r = 0.0, ti = dat_scramble$ti, ar1 = 0.5)
  V1s <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, smooth_vi = TRUE)
  V2s <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti, smooth_vi = TRUE)
  V3s <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.6, ti = dat_scramble$ti, ar1 = 0, smooth_vi = TRUE)
  V4s <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, ti = dat_scramble$ti, ar1 = 0.5, smooth_vi = TRUE)
  V5s <- impute_covariance_matrix(vi = dat_scramble$vi, cluster = dat_scramble$study, r = 0.0, ti = dat_scramble$ti, ar1 = 0.5, smooth_vi = TRUE)
  
  expect_equal(V1, V1s)
  expect_equal(V2, V2s)
  expect_equal(V3, V3s)
  expect_equal(V4, V4s)
  expect_equal(V5, V5s)
  
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
  expect_equal(V_resorted, unblock(V_unscramble))
})

test_that("impute_covariance_matrix works with missing variances.", {
  skip_if_not_installed("robumeta")
  
  data(corrdat, package = "robumeta")
  dat_miss <- corrdat
  dat_miss$var[sample.int(nrow(corrdat), size = round(nrow(corrdat) / 10))] <- NA
  V_missing <- impute_covariance_matrix(dat_miss$var, cluster = dat_miss$studyid, 
                                        r = 0.8, return_list = FALSE)
  non_missing_rows <- !is.na(dat_miss$var)
  V_missing <- V_missing[non_missing_rows, non_missing_rows]
  
  dat_complete <- subset(dat_miss, !is.na(var))
  V_complete <- impute_covariance_matrix(dat_complete$var, cluster = dat_complete$studyid, 
                                         r = 0.8, return_list = FALSE)
  
  expect_equal(V_missing, V_complete)
})




test_that("pattern_covariance_matrix works.", {
  skip_if_not_installed("robumeta")
  skip_if_not_installed("metafor")
  
  data(oswald2013, package = "robumeta")
  dat <- metafor::escalc(data = oswald2013, measure = "ZCOR", ri = R, ni = N)
  
  # make a patterned correlation matrix 
  p_levels <- levels(dat$Crit.Cat)
  r_pattern <- 0.7^as.matrix(dist(1:length(p_levels)))
  diag(r_pattern) <- seq(0.75, 0.95, length.out = 6)
  rownames(r_pattern) <- colnames(r_pattern) <- p_levels
  
  # impute the covariance matrix using patterned correlations
  V_list <- pattern_covariance_matrix(vi = dat$vi, 
                                      cluster = dat$Study, 
                                      pattern_level = dat$Crit.Cat,
                                      r_pattern = r_pattern,
                                      smooth_vi = TRUE)

  r_pattern1 <- r_pattern2 <- r_pattern 
  
  # Recreate a constant covariance matrix
  
  r_pattern1[1:6,1:6] <- 0.6
  
  V_pattern_const <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                              r_pattern = r_pattern1, smooth_vi = FALSE)
  
  V_impute_const <- impute_covariance_matrix(vi = dat$vi, cluster = dat$Study, r = 0.6)
  expect_equal(V_pattern_const, V_impute_const, check.attributes = FALSE)
  
  # Patterns work with excluded categories
  exclude <- 3:5
  r_pattern2[3:5,] <- 0.3
  r_pattern2[,3:5] <- 0.3
  
  V_pattern_exclude <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                                 r_pattern = r_pattern[-(3:5),-(3:5)], r = 0.3, smooth_vi = FALSE)
  V_pattern_full <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                                 r_pattern = r_pattern2, smooth_vi = FALSE)
  expect_equal(V_pattern_exclude, V_pattern_full)
  
  # Patterns work with extra categories
  levs <- c(p_levels, LETTERS[1:5])
  r_pattern3 <- matrix(-4, nrow = length(levs), ncol = length(levs), dimnames = list(levs, levs))
  r_pattern3[p_levels,p_levels] <- r_pattern[p_levels, p_levels]
  
  V_pattern_extra <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                                 r_pattern = r_pattern3, r = 0.3, smooth_vi = TRUE)
  expect_equal(V_pattern_extra, V_list)
  
  V_pattern_extra_exclude <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                                      r_pattern = r_pattern3[-(3:5),-(3:5)], r = 0.3, smooth_vi = FALSE)
  expect_equal(V_pattern_extra_exclude, V_pattern_full)
  
  VS_pattern_extra_exclude <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                                       r_pattern = r_pattern3[-(3:5),-(3:5)], r = 0.3, 
                                                       subgroup = dat$IAT.Focus,
                                                       smooth_vi = FALSE, return_list = FALSE)
  VS_pattern_full <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                              r_pattern = r_pattern2, 
                                              subgroup = dat$IAT.Focus,
                                              smooth_vi = FALSE, return_list = FALSE)
  
  expect_equal(VS_pattern_extra_exclude, VS_pattern_full)
  
  
  expect_warning(
    pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                              r_pattern = r_pattern3[-(3:5),-(3:5)], r = 4, 
                              subgroup = dat$IAT.Focus,
                              smooth_vi = FALSE, return_list = FALSE)
  )

  V_npd <- pattern_covariance_matrix(vi = dat$vi, cluster = dat$Study, pattern_level = dat$Crit.Cat,
                                     r_pattern = r_pattern3[-(3:5),-(3:5)], r = 4, 
                                     subgroup = dat$IAT.Focus,
                                     smooth_vi = FALSE, return_list = FALSE, check = FALSE)  
  expect_true(inherits(V_npd, "matrix"))

  # pattern_covariance_matrix works with missing entries
  dat_miss <- dat
  dat_miss$vi[sample.int(nrow(dat), size = round(nrow(dat) / 10))] <- NA
  
  V_missing <- pattern_covariance_matrix(dat_miss$vi, cluster = dat_miss$Study,
                                         pattern_level = dat_miss$Crit.Cat,
                                         r_pattern = r_pattern, return_list = FALSE)
  
  non_missing_rows <- !is.na(dat_miss$vi)
  V_missing <- V_missing[non_missing_rows, non_missing_rows]
  
  dat_complete <- subset(dat_miss, !is.na(vi))
  V_complete <- pattern_covariance_matrix(dat_complete$vi, cluster = dat_complete$Study, 
                                          pattern_level = dat_complete$Crit.Cat,
                                          r_pattern = r_pattern, return_list = FALSE)
  
  expect_equal(V_missing, V_complete)
  
  
  dat_miss$Crit.Cat[sample.int(nrow(dat), size = round(nrow(dat) / 10))] <- NA
  
  expect_error(
    pattern_covariance_matrix(dat_miss$vi, cluster = dat_miss$Study,
                              pattern_level = dat_miss$Crit.Cat,
                              r_pattern = r_pattern, return_list = TRUE,
                              check_PD = FALSE)
  )  
})

