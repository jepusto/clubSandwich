context("lm objects")

m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
n <- length(cluster)
X <- matrix(rnorm(3 * n), n, 3)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
w <- rgamma(n, shape = 3, scale = 3)
y <- X %*% c(.4, .3, -.3) + nu + e

dat <- data.frame(y, X, cluster, w, row = 1:n)

lm_fit <- lm(y ~ X1 + X2 + X3, data = dat)
WLS_fit <- lm(y ~ X1 + X2 + X3, data = dat, weights = w)
CR_types <- paste0("CR",0:4)

# obj <- WLS_fit
# type <- "CR2" 
# vcov <- vcovCR(obj, cluster = cluster, type = type)
# target = NULL
# inverse_var = FALSE

test_that("vcovCR options don't matter for CR0", {
  expect_error(vcovCR(lm_fit, type = "CR0"))
  CR0 <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0")
  attr(CR0, "target") <- NULL
  attr(CR0, "inverse_var") <- NULL
  CR0_A <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w)
  attr(CR0_A, "target") <- NULL
  attr(CR0_A, "inverse_var") <- NULL
  expect_identical(CR0_A, CR0)
  CR0_B <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = FALSE)
  attr(CR0_B, "target") <- NULL
  attr(CR0_B, "inverse_var") <- NULL
  expect_identical(CR0_A, CR0)
  CR0_C <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = TRUE)
  attr(CR0_C, "target") <- NULL
  attr(CR0_C, "inverse_var") <- NULL
  expect_identical(CR0_C, CR0)
  
  wCR0 <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0")
  attr(wCR0, "target") <- NULL
  attr(wCR0, "inverse_var") <- NULL
  wCR0_A <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w)
  attr(wCR0_A, "target") <- NULL
  attr(wCR0_A, "inverse_var") <- NULL
  expect_identical(wCR0_A, wCR0)
  wCR0_B <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = FALSE)
  attr(wCR0_B, "target") <- NULL
  attr(wCR0_B, "inverse_var") <- NULL
  expect_identical(wCR0_B, wCR0)
  wCR0_C <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = TRUE)
  attr(wCR0_C, "target") <- NULL
  attr(wCR0_C, "inverse_var") <- NULL
  expect_identical(wCR0_C, wCR0)
})

test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)

  attr(CR2_iv, "inverse_var") <- FALSE
  CR2_not <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = FALSE)
  expect_equal(CR2_not, CR2_iv)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n)), CR2_not)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w), CR2_not))
  
  wCR2_id <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", inverse_var = FALSE), wCR2_id)
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n)), wCR2_id)
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = FALSE), wCR2_id)
  
  wCR2_iv <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
  wCR2_target <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE)
  expect_false(identical(wCR2_target, wCR2_id))
  expect_identical(matrix(wCR2_target, dim(wCR2_target)), matrix(wCR2_iv, dim(wCR2_iv)))
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE), wCR2_target)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR4")
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = TRUE), CR4_iv)
  
  attr(CR4_iv, "inverse_var") <- FALSE
  CR4_not <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", inverse_var = FALSE)
  expect_equal(CR4_not, CR4_iv)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n)), CR4_not)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = 1 / dat$w), CR4_not))
  
  wCR4_id <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4")
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", inverse_var = FALSE), wCR4_id)
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n)), wCR4_id)
  expect_identical(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = FALSE), wCR4_id)
  
  wCR4_iv <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", inverse_var = TRUE)
  wCR4_target <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = 1 / dat$w, inverse_var = TRUE)
  expect_false(identical(wCR4_target, wCR4_id))
  expect_identical(matrix(wCR4_target, dim(wCR4_target)), matrix(wCR4_iv, dim(wCR4_iv)))
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = 1 / dat$w, inverse_var = TRUE), wCR4_target)
  
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(lm_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(WLS_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(lm_fit, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(WLS_fit, vcov = "CR4", cluster = dat$cluster))
})

test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)
  CR0 <- vcovCR(lm_fit, cluster = dat$row, type = "CR0")
  expect_equal(vcovHC(lm_fit, type = "HC0"), as.matrix(CR0))
  CR1 <- vcovCR(lm_fit, cluster = dat$row, type = "CR1S")
  expect_equal(vcovHC(lm_fit, type = "HC1"), as.matrix(CR1))
  CR2 <- vcovCR(lm_fit, cluster = dat$row, type = "CR2")
  expect_equal(vcovHC(lm_fit, type = "HC2"), as.matrix(CR2))
  CR3 <- vcovCR(lm_fit, cluster = dat$row, type = "CR3")
  expect_equal(vcovHC(lm_fit, type = "HC3"), as.matrix(CR3))
})

test_that("CR2 is equivalent to Welch t-test for DiD design", {
  m0 <- 4
  m1 <- 9
  m <- m0 + m1
  cluster <- factor(rep(LETTERS[1:m], each = 2))
  n <- length(cluster)
  time <- rep(c(1,2), m)
  trt_clusters <- c(rep(0,m0), rep(1,m1))
  trt <- (time - 1) * rep(trt_clusters, each = 2)
  nu <- rnorm(m)[cluster]
  e <- rnorm(n)
  y <- 0.4 * trt + nu + e
  
  dat <- data.frame(y, time, trt, cluster)
  lm_DID <- lm(y ~ cluster + factor(time) + trt, data = dat)
  t_Satt <- coef_test(lm_DID, vcov = "CR2", cluster = dat$cluster)["trt",]
  y_diff <- apply(matrix(y, nrow = 2), 2, diff)
  t_Welch <- t.test(y_diff ~ trt_clusters)
  
  expect_equal(with(t_Welch, estimate[[2]] - estimate[[1]]), t_Satt$beta)
  expect_equal(as.numeric(-t_Welch$statistic), with(t_Satt, beta / SE))
  expect_is(all.equal(as.numeric(t_Welch$parameter), t_Satt$df), "character")
})

test_that("Order doesn't matter.",{
  dat_scramble <- dat[sample(n),]
  WLS_scramble <- update(WLS_fit, data = dat_scramble)

  CR_fit <- lapply(CR_types, function(x) vcovCR(WLS_fit, cluster = dat$cluster, type = x))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(WLS_scramble, cluster = dat_scramble$cluster, type = x))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(WLS_fit, vcov = x, cluster = dat$cluster, test = "All"))
  test_scramble <- lapply(CR_types, function(x) coef_test(WLS_scramble, vcov = x, cluster = dat_scramble$cluster, test = "All"))
  expect_equal(test_fit, test_scramble, tolerance = 10^-6)
  
  constraints <- combn(length(coef(lm_fit)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(WLS_fit, constraints = constraints, vcov = "CR2", cluster = dat$cluster, test = "All")
  Wald_scramble <- Wald_test(WLS_scramble, constraints = constraints, vcov = "CR2", cluster = dat_scramble$cluster, test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})

test_that("clubSandwich works with dropped observations", {
  dat_miss <- dat
  dat_miss$X1[sample.int(n, size = round(n / 10))] <- NA
  lm_dropped <- lm(y ~ X1 + X2 + X3, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(X1))
  lm_complete <- lm(y ~ X1 + X2 + X3, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, cluster = dat_miss$cluster, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, cluster = dat_complete$cluster, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, cluster = dat_miss$cluster, test = "All"))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, cluster = dat_complete$cluster, test = "All"))
  expect_identical(test_drop, test_complete)
})
