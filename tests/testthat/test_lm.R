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

lm_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat)
WLS_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat, weights = w)

obj <- WLS_fit
type <- "CR2" 
target = NULL
inverse_var = FALSE

test_that("vcovCR options don't matter for CR0", {
  expect_error(vcovCR(lm_fit, type = "CR0"))
  CR0 <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0")
  attr(CR0, "target") <- NULL
  CR0_A <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w)
  attr(CR0_A, "target") <- NULL
  expect_identical(CR0_A, CR0)
  CR0_B <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = FALSE)
  attr(CR0_B, "target") <- NULL
  expect_identical(CR0_A, CR0)
  CR0_C <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = TRUE)
  attr(CR0_C, "target") <- NULL
  expect_identical(CR0_C, CR0)
  
  wCR0 <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0")
  attr(wCR0, "target") <- NULL
  wCR0_A <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w)
  attr(wCR0_A, "target") <- NULL
  expect_identical(wCR0_A, wCR0)
  wCR0_B <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = FALSE)
  attr(wCR0_B, "target") <- NULL
  expect_identical(wCR0_B, wCR0)
  wCR0_C <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0", target = 1 / dat$w, inverse_var = TRUE)
  attr(wCR0_C, "target") <- NULL
  expect_identical(wCR0_C, wCR0)
})

test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)

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
#  expect_identical(matrix(wCR2_target, dim(wCR2_target)), matrix(wCR2_iv, dim(wCR2_iv)))
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE), wCR2_target)
  
})

test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich)
  CR0 <- vcovCR(lm_fit, cluster = dat$row, type = "CR0")
  expect_equal(vcovHC(lm_fit, type = "HC0"), as.matrix(CR0))
  CR1 <- vcovCR(lm_fit, cluster = dat$row, type = "CR1S")
  expect_equal(vcovHC(lm_fit, type = "HC1"), as.matrix(CR1) * (n - 1) / n)
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
