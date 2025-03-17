context("confidence intervals")
set.seed(20190513)

skip_if_not_installed("nlme")

library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Ovary, package = "nlme")

Ovary$time_int <- 1:nrow(Ovary)

gls_fit <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
              correlation = corAR1(form = ~ time_int | Mare), 
              weights = varPower())

CRs <- paste0("CR", 0:4)

test_that("vcov arguments work", {
  VCR <- lapply(CRs, function(t) vcovCR(gls_fit, type = t))
  CI_A <- lapply(VCR, function(v) conf_int(gls_fit, vcov = v, level = .98))
  CI_B <- lapply(CRs, function(t) conf_int(gls_fit, vcov = t, level = .98))
  expect_equal(CI_A, CI_B)
})

test_that("coefs argument works", {
  which_grid <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(gls_fit))))
  tests_all <- conf_int(gls_fit, vcov = "CR0", coefs = "All")
  
  CI_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
  CI_B <- apply(which_grid[-1,], 1, function(x) conf_int(gls_fit, vcov = "CR0", coefs = x))
  expect_equal(CI_A, CI_B)
})

test_that("printing works", {
  CIs <- conf_int(gls_fit, vcov = "CR0")
  expect_output(print(CIs))
  
  CIs <- conf_int(gls_fit, vcov = "CR0", p_values = TRUE)
  expect_output(x <- print(CIs))
  expect_true(all(c("p-value","Sig.") %in% names(x)))
})

test_that("level checks work", {
  expect_error(conf_int(gls_fit, vcov = "CR0", level = -0.01))
  expect_error(conf_int(gls_fit, vcov = "CR0", level = 95))
  expect_output(print(conf_int(gls_fit, vcov = "CR0", level = runif(1))))
})

test_that("CI boundaries are ordered", {
  lev <- runif(1)
  CI_z <- conf_int(gls_fit, vcov = "CR0", test = "z", level = lev)
  CI_t <- conf_int(gls_fit, vcov = "CR0", test = "naive-t", level = lev)
  CI_Satt <- conf_int(gls_fit, vcov = "CR0", test = "Satterthwaite", level = lev)
  expect_true(all(CI_t$CI_L < CI_z$CI_L))
  expect_true(all(CI_t$CI_U > CI_z$CI_U))
  expect_true(all(CI_Satt$CI_L < CI_z$CI_L))
  expect_true(all(CI_Satt$CI_U > CI_z$CI_U))
})

test_that("conf_int() is consistent with coef_test()", {
  
  lev <- runif(1)
  CIs <- lapply(CRs, function(v) conf_int(gls_fit, vcov = v, test = "Satterthwaite", level = lev, p_values = TRUE))
  ttests <- lapply(CRs, function(v) coef_test(gls_fit, vcov = v, test = "Satterthwaite"))
  CI_L <- lapply(ttests, function(x) x$beta - x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  CI_U <- lapply(ttests, function(x) x$beta + x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  expect_equal(lapply(CIs, function(x) x$CI_L), CI_L)
  expect_equal(lapply(CIs, function(x) x$CI_U), CI_U)
  expect_equal(lapply(CIs, function(x) x$p_val), lapply(ttests, function(x) x$p_Satt))

  lev <- runif(1)
  CIs <- lapply(CRs, function(v) conf_int(gls_fit, vcov = v, test = "naive-t", level = lev, p_values = TRUE))
  ttests <- lapply(CRs, function(v) coef_test(gls_fit, vcov = v, test = "naive-t"))
  CI_L <- lapply(ttests, function(x) x$beta - x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  CI_U <- lapply(ttests, function(x) x$beta + x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  expect_equal(lapply(CIs, function(x) x$CI_L), CI_L)
  expect_equal(lapply(CIs, function(x) x$CI_U), CI_U)
  expect_equal(lapply(CIs, function(x) x$p_val), lapply(ttests, function(x) x$p_t))
  
  lev <- runif(1)
  CIs <- lapply(CRs, function(v) conf_int(gls_fit, vcov = v, test = "z", level = lev, p_values = TRUE))
  ttests <- lapply(CRs, function(v) coef_test(gls_fit, vcov = v, test = "z"))
  CI_L <- lapply(ttests, function(x) x$beta - x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  CI_U <- lapply(ttests, function(x) x$beta + x$SE * qt(1 - (1 - lev) / 2, df = x$df))
  expect_equal(lapply(CIs, function(x) x$CI_L), CI_L)
  expect_equal(lapply(CIs, function(x) x$CI_U), CI_U)
  expect_equal(lapply(CIs, function(x) x$p_val), lapply(ttests, function(x) x$p_z))
  
})

test_that("conf_int has informative error messages.", {
  expect_error(
    conf_int(gls_fit, vcov = "CR0", test = "all")
  )
  
  expect_error(
    conf_int(gls_fit, vcov = "CR0", test = "saddlepoint")
  )
})

test_that("linear_contrast multiple comparisons p-value adjustment works correctly", {
  
  # taken from example usage
  data("ChickWeight", package = "datasets")
  lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
  
  lc1 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick, 
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE)
  
  lc2 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE, adjustment_method = "none")

  expect_equal(lc1, lc2) # check explicitly stating adjustment_method default doesn't change anything
  
  # test using adjustment without p_values = TRUE
  expect_warning(lc1 <- linear_contrast(lm_fit,
                                        vcov = "CR2",
                                        cluster = ChickWeight$Chick,
                                        contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                                        p_values = FALSE, adjustment_method = "hochberg"))
  
  # check that the above code performs default behavior
  lc2 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = FALSE)
  expect_equal(lc1, lc2)

  # lc2 will now show p_values with no adjustment
  lc2 <- linear_contrast(lm_fit,
                         vcov = "CR2",
                         cluster = ChickWeight$Chick,
                         contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                         p_values = TRUE)
  
  # test using nonexistent adjustment type
  expect_warning(lc1 <- linear_contrast(lm_fit,
                                        vcov = "CR2",
                                        cluster = ChickWeight$Chick, 
                                        contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
                                        p_values = TRUE,
                                        adjustment_method = "nonexistent"),
                 "The specified adjustment method")
  # check that the above defaults to no adjustment
  expect_equal(lc1, lc2)
  
  # test using p-adjustment when results have a length of 1
  # TBD, haven't figured out how to get something with only one row
  
})