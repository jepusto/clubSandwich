context("t-tests")

balanced_dat <- function(m, n) {
  cluster <- factor(rep(LETTERS[1:m], each = n))
  N <- length(cluster)
  X_btw <- rep(rep(LETTERS[1:3], c(3,5,7)), each = n)
  X_wth <- rep(rep(c(0,1), each = n / 2), m)
  nu <- rnorm(m)[cluster]
  e <- rnorm(n)
  y <- nu + e
  data.frame(y, X_btw, X_wth, cluster, row = 1:N)
  
}

CRs <- paste0("CR", 0:4)

test_that("vcov arguments work", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  VCR <- lapply(CRs, function(t) vcovCR(lm_fit, cluster = dat$cluster, type = t))
  test_A <- lapply(VCR, function(v) coef_test(lm_fit, vcov = v, test = "All"))
  test_B <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, cluster = dat$cluster, test = "All"))
  expect_identical(test_A, test_B)
})

test_that("printing works", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  test_results <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, cluster = dat$cluster, test = "All"))
  expect_identical(test_A, test_B)
})

test_that("p-values are ordered", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  test_results <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, cluster = dat$cluster, test = "All"))
  expect_identical(test_A, test_B)
})
