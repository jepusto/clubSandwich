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
  expect_equal(test_A, test_B)
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

test_that("coefs argument works", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  which_grid <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(lm_fit))))
  tests_all <- coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
  
  tests_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
  tests_B <- apply(which_grid[-1,], 1, function(x) coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = x, p_values = FALSE))
  expect_equal(tests_A, tests_B)
})

test_that("printing works", {
  dat <- balanced_dat(m = 15, n = 8)
  lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
  t_tests <- coef_test(lm_fit, vcov = "CR2", cluster = dat$cluster, test = "All")
  expect_output(print(t_tests))
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
