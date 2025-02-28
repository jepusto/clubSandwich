context("t-tests")
set.seed(20190513)

balanced_dat <- function(m, n) {
  cluster <- factor(rep(LETTERS[1:m], each = n))
  N <- length(cluster)
  m1 <- sample(3:(m-7), size = 1)
  m2 <- sample((m1 + 3):(m-3), size = 1) - m1
  m3 <- m - m1 - m2
  c(m1, m2, m3)
  X_btw <- rep(rep(LETTERS[1:3], c(m1, m2, m3)), each = n)
  X_wth <- rep(rep(c(0,1), each = n / 2), m)
  nu <- rnorm(m)[cluster]
  e <- rnorm(n * m)
  y <- nu + e
  data.frame(y, X_btw, X_wth, cluster, row = 1:N)
}

CRs <- paste0("CR", 0:4)

test_that("vcov arguments work", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  VCR <- lapply(CRs, function(t) vcovCR(lm_fit, cluster = dat$cluster, type = t))
  test_A <- lapply(VCR, function(v) coef_test(lm_fit, vcov = v, test = "All", p_values = FALSE))
  test_B <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, cluster = dat$cluster, test = "All", p_values = FALSE))
  test_C <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, null_constant = 0, cluster = dat$cluster, test = "All", p_values = FALSE))
  compare_ttests(test_A, test_B)
  compare_ttests(test_A, test_C)
})

test_that("get_which_coef() works", {
  f <- 6
  beta <- 1:f
  beta_names <- letters[1:f]
  names(beta) <- beta_names
  
  which_grid <- as.matrix(expand.grid(rep(list(c(FALSE,TRUE)), f)))
  dimnames(which_grid) <- NULL
  name_list <- apply(which_grid, 1, function(x) beta_names[x])
  int_list <- apply(which_grid, 1, which)
  
  which_log <- apply(which_grid[-1,], 1, get_which_coef, beta = beta)
  which_char <- sapply(name_list[-1], get_which_coef, beta = beta)
  which_int <- sapply(int_list[-1], get_which_coef, beta = beta)
  
  expect_identical(get_which_coef(beta, coefs = "All"), rep(TRUE, f))
  expect_error(get_which_coef(beta, coefs = which_grid[1,]))
  expect_error(get_which_coef(beta, coefs = name_list[[1]]))
  expect_error(get_which_coef(beta, coefs = int_list[[1]]))
  
  expect_identical(which_log, which_char)
  expect_identical(which_log, which_int)
  expect_identical(which_char, which_int)
})

test_that("coefs argument works.", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  which_grid <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(lm_fit))))
  tests_all <- coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
  
  tests_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
  tests_B <- apply(which_grid[-1,], 1, function(x) coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = x, p_values = FALSE))
  expect_equal(tests_A, tests_B, check.attributes = FALSE)
})

test_that("coefs argument works with null_constants.", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  which_grid <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(lm_fit))))

  tests_all <- coef_test(lm_fit, vcov = "CR0", null_constants = 0.1, cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
  tests_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
  tests_B <- apply(which_grid[-1,], 1, function(x) coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = x, null_constants = 0.1, p_values = FALSE))
  expect_equal(tests_A, tests_B, check.attributes = FALSE)

  null_consts <- rnorm(4)
  tests_all <- coef_test(lm_fit, vcov = "CR0", null_constants = null_consts, cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
  tests_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
  tests_B <- apply(which_grid[-1,], 1, function(x) coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = x, null_constants = null_consts[x], p_values = FALSE))
  expect_equal(tests_A, tests_B, check.attributes = FALSE)
  
})

test_that("printing works", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  t_tests <- coef_test(lm_fit, vcov = "CR2", cluster = dat$cluster, test = "All")
  expect_output(x <- print(t_tests))

  expect_equal(t_tests$df_z, rep(Inf, 4L))
  expect_equal(t_tests$df_t, rep(14L, 4L))
  expect_true(all(t_tests$df_t >= round(t_tests$df_Satt,1)))
  
  expect_identical(names(x),
                   c("Coef.","Estimate","SE","t-stat",
                     "d.f. (z)", "p-val (z)", "Sig.",
                     "d.f. (naive-t)", "p-val (naive-t)","Sig.",
                     "d.f. (naive-tp)", "p-val (naive-tp)","Sig.",
                     "d.f. (Satt)", "p-val (Satt)", "Sig.",
                     "s.p.", "p-val (Saddle)", "Sig."))
})

test_that("p-values are ordered", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  test_results <- lapply(CRs, function(t) coef_test(lm_fit, vcov = t, cluster = dat$cluster, test = "All"))
  test_results <- do.call(rbind, test_results)
  expect_true(with(test_results, all(p_z < p_t)))
  expect_true(with(test_results, all(p_z < p_Satt)))
})

test_that("Satterthwaite df work for special cases", {
  m <- sample(12:26, size = 1)
  n <- sample(seq(4,12,2), size = 1)
  dat <- balanced_dat(m, n)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  t_tests <- coef_test(lm_fit, vcov = "CR2", cluster = dat$cluster, test = "Satterthwaite")
  expect_equal(t_tests$df[4], m - 1)
  mg <- table(dat$X_btw) / n
  df <- apply(cbind(mg[1], mg[-1]), 1, function(x) sum(x)^2 * prod(x - 1) / sum(x^2 * (x - 1)))
  expect_equivalent(t_tests$df[2:3], df)
  lm_fit <- lm(y ~ 0 + cluster + X_wth, data = dat)
  t_tests <- coef_test(lm_fit, vcov = "CR2", cluster = dat$cluster, test = "Satterthwaite")
  expect_equal(t_tests$df[m + 1], m - 1) 
})

test_that("null_constants error messages appear as appropriate.", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  vcr <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0")
  expect_error(
    coef_test(lm_fit, vcov = vcr, null_constants = 1:3)
  )
  expect_error(
    coef_test(lm_fit, vcov = vcr, coefs = 1:2, null_constants = rep(0,3))
  )
  expect_error(
    coef_test(lm_fit, vcov = vcr, null_constants = "none")
  )
})
