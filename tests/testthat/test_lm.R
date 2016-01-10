context("lm objects")

m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
table(cluster)
n <- length(cluster)
X <- matrix(rnorm(3 * n), n, 3)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
w <- rgamma(n, shape = 3, scale = 3)
y <- X %*% c(.4, .3, -.3) + nu + e

dat <- data.frame(y, X, cluster, w)

lm_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat)
WLS_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat, weights = w)

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
  expect_identical(matrix(wCR2_target, dim(wCR2_target)), matrix(wCR2_iv, dim(wCR2_iv)))
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE), wCR2_target)
  
})

# test for equality with HC when cluster = rownames(data)
