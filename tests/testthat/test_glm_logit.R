context("logit glm objects")

m <- 20
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
n <- length(cluster)
X1 <- c(rep(-0.5, m / 2), rep(0.5, m / 2))[cluster]
X2 <- c(rep(-0.3, 0.4 * m), rep(0.7, 0.3 * m), rep(-0.3, 0.4 * m))[cluster]
X3 <- rnorm(m)[cluster] + rnorm(n)
X4 <- rnorm(n)
X <- cbind(X1, X2, X3, X4)
eta <- -0.4 + X %*% c(0.3, -0.6, 0.15, 0.15)
p <- 1 / (1 + exp(-eta))
summary(p)

w <- sample(1:4, size = n, replace = TRUE)
y1 <- rbinom(n, size = 1, prob = p)
y2 <- rbinom(n, size = w, prob = p)
yp <- y2 / w

dat <- data.frame(y1, y2, yp, X, cluster, w, row = 1:n)

logit_fit <- glm(y1 ~ X1 + X2 + X3 + X4, data = dat, family = "binomial")
sflogit_fit <- glm(cbind(y2, w - y2) ~ X1 + X2 + X3 + X4, data = dat, family = "binomial")
plogit_fit <- glm(yp ~ X1 + X2 + X3 + X4, data = dat, weights = w, family = "quasibinomial")

# obj <- logit_fit
# y <- dat$y1
# type <- "CR2"
# vcov <- vcovCR(obj, cluster = cluster, type = type)
# target = NULL
# inverse_var = FALSE
# 
# cluster <- droplevels(as.factor(cluster))
# B <- sandwich::bread(obj) / v_scale(obj)
# X_list <- matrix_list(model_matrix(obj), cluster, "row")
# W_list <- weightMatrix(obj, cluster)
# XWX <- Reduce("+", Map(function(x, w) t(x) %*% w %*% x, x = X_list, w = W_list))
# M <- chol2inv(chol(XWX))
# attr(M, "dimnames") <- attr(B, "dimnames")
# 
# M / B
# diff(range(M / B))


test_that("bread works", {
  
  expect_true(check_bread(logit_fit, cluster = dat$cluster, check_coef = FALSE, tol = 10^-3))
  glm_vcov <- bread(logit_fit) * summary(logit_fit)$dispersion / v_scale(logit_fit)
  expect_equal(vcov(logit_fit), glm_vcov)
  
  expect_true(check_bread(sflogit_fit, cluster = dat$cluster, check_coef = FALSE, tol = 10^-3))
  glm_vcov <- bread(sflogit_fit) * summary(sflogit_fit)$dispersion / v_scale(sflogit_fit)
  expect_equal(vcov(sflogit_fit), glm_vcov)

  expect_true(check_bread(plogit_fit, cluster = dat$cluster, check_coef = FALSE, tol = 10^-3))
  glm_vcov <- bread(plogit_fit) * summary(plogit_fit)$dispersion / v_scale(plogit_fit)
  expect_equal(vcov(plogit_fit), glm_vcov)
  
})

test_that("vcovCR options work for CR2", {
  
  CR2_iv <- vcovCR(logit_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(logit_fit, cluster = dat$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR2", 
                          target = targetVariance(logit_fit, cluster = dat$cluster), 
                          inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)

  CR2_iv <- vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)
  
  CR2_iv <- vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2")
  expect_identical(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(logit_fit, cluster = dat$cluster, type = "CR4")
  expect_identical(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
  CR4_iv <- vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4")
  expect_identical(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
  CR4_iv <- vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4")
  expect_identical(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
})

test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(logit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(sflogit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(plogit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(logit_fit, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(sflogit_fit, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(plogit_fit, vcov = "CR4", cluster = dat$cluster))
})

test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)
  
  HC_types <- paste0("HC", 0:3)
  HC_list <- lapply(HC_types, function(t) vcovHC(logit_fit, type = t))
  
  CR_types <- paste0("CR", 0:3)
  CR_types[2] <- "CR1S"
  CR_list <- lapply(CR_types, function(t) as.matrix(vcovCR(logit_fit, cluster = dat$row, type = t)))
  expect_equal(HC_list, CR_list, tol = 4 * 10^-4)
  
})

CR_types <- paste0("CR", 0:4)

test_that("Order doesn't matter.",{
  dat_scramble <- dat[sample(n),]
  logit_scramble <- update(logit_fit, data = dat_scramble)

  CR_fit <- lapply(CR_types, function(x) vcovCR(logit_fit, cluster = dat$cluster, type = x))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(logit_scramble, cluster = dat_scramble$cluster, type = x))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(logit_fit, vcov = x, cluster = dat$cluster, test = "All"))
  test_scramble <- lapply(CR_types, function(x) coef_test(logit_scramble, vcov = x, cluster = dat_scramble$cluster, test = "All"))
  compare_tests <- mapply(function(a, b) max(abs(a / b - 1), na.rm = TRUE), test_fit, test_scramble)
  expect_true(all(compare_tests < 10^-6))
  
  constraints <- combn(length(coef(logit_fit)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(logit_fit, constraints = constraints, vcov = "CR2", cluster = dat$cluster, test = "All")
  Wald_scramble <- Wald_test(logit_scramble, constraints = constraints, vcov = "CR2", cluster = dat_scramble$cluster, test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})


test_that("clubSandwich works with dropped observations", {
  dat_miss <- dat
  dat_miss$X1[sample.int(n, size = round(n / 10))] <- NA
  logit_dropped <- update(logit_fit, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(X1))
  logit_complete <- update(logit_fit, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(logit_dropped, cluster = dat_miss$cluster, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(logit_complete, cluster = dat_complete$cluster, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(logit_dropped, vcov = x, cluster = dat_miss$cluster, test = "All"))
  test_complete <- lapply(CR_types, function(x) coef_test(logit_complete, vcov = x, cluster = dat_complete$cluster, test = "All"))
  expect_identical(test_drop, test_complete)
})


test_that("clubSandwich works with aliased predictors", {
  data(npk, package = "datasets")
  npk_alias <- glm(yield ~ block + N*P*K, data = npk)
  npk_drop <- glm(yield ~ block + N + P + K + N:P + N:K + P:K, data = npk)
  
  CR_alias <- lapply(CR_types[-4], function(x) vcovCR(npk_alias, cluster = npk$block, type = x))
  CR_drop <- lapply(CR_types[-4], function(x) vcovCR(npk_drop, cluster = npk$block, type = x))
  expect_identical(CR_alias, CR_drop)
  
  test_drop <- lapply(CR_types[-4], function(x) coef_test(npk_alias, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite"))[-13,])
  test_complete <- lapply(CR_types[-4], function(x) coef_test(npk_drop, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite")))
  expect_identical(test_drop, test_complete)
})


test_that("clubSandwich results are equivalent to geepack", {
  library(geepack)
  
  # check CR0 with logit
  logit_gee <- geeglm(y1 ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, family = "binomial")
  logit_refit <- update(logit_fit, start = coef(logit_gee))
  expect_equal(coef(logit_refit), coef(logit_gee))
  
  V_gee0 <- summary(logit_gee)$cov.scaled
  V_CR0 <- as.matrix(vcovCR(logit_refit, cluster = dat$cluster, type = "CR0"))
  attr(V_gee0, "dimnames") <- attr(V_CR0, "dimnames")
  expect_equal(V_gee0, V_CR0)
  
  # check CR3 with logit
  logit_gee <- geeglm(y1 ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, family = "binomial",
                      std.err = "jack")
  V_gee3 <- summary(logit_gee)$cov.scaled
  V_CR3 <- as.matrix(vcovCR(logit_refit, cluster = dat$cluster, type = "CR3"))
  attr(V_gee3, "dimnames") <- attr(V_CR3, "dimnames")
  expect_equal(V_gee3 * m / (m - 6), V_CR3)
  
  # check CR0 with plogit
  plogit_gee <- geeglm(yp ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, weights = w, family = "binomial")
  plogit_refit <- update(plogit_fit, start = coef(plogit_gee))
  expect_equal(coef(plogit_refit), coef(plogit_gee))
  
  V_gee0 <- summary(plogit_gee)$cov.scaled
  V_CR0 <- as.matrix(vcovCR(plogit_refit, cluster = dat$cluster, type = "CR0"))
  attr(V_gee0, "dimnames") <- attr(V_CR0, "dimnames")
  expect_equal(V_gee0, V_CR0)
  
})
