context("lm objects")
set.seed(20190513)

m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
n <- length(cluster)
const <- rep("whuzzup.", n)
X <- matrix(rnorm(3 * n), n, 3)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
w <- rgamma(n, shape = 3, scale = 3)
y <- X %*% c(.4, .3, -.3) + nu + e

dat <- data.frame(y, X, cluster, const, w, row = 1:n)

lm_fit <- lm(y ~ X1 + X2 + X3, data = dat)
WLS_fit <- lm(y ~ X1 + X2 + X3, data = dat, weights = w)
CR_types <- paste0("CR",0:4)

# obj <- WLS_fit
# y <- dat$y
# type <- "CR2"
# vcov <- vcovCR(obj, cluster = cluster, type = type)
# target = NULL
# inverse_var = FALSE

test_that("bread works", {
  
  expect_true(check_bread(lm_fit, cluster = dat$cluster, y = dat$y))
  lm_vcov <- bread(lm_fit) * summary(lm_fit)$sigma^2 / v_scale(lm_fit)
  expect_equal(vcov(lm_fit), lm_vcov)
  
  expect_true(check_bread(WLS_fit, cluster = dat$cluster, y = dat$y))
  wls_vcov <- bread(WLS_fit) * summary(WLS_fit)$sigma^2 / v_scale(WLS_fit)
  expect_equal(vcov(WLS_fit), wls_vcov)
})

test_that("vcovCR options don't matter for CR0", {
  expect_error(vcovCR(lm_fit, type = "CR0"))
  CR0 <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0")
  expect_output(print(CR0))
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
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)
  
  attr(CR2_iv, "inverse_var") <- FALSE
  CR2_not <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = FALSE)
  expect_equal(CR2_not, CR2_iv)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n)), CR2_not)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w), CR2_not))
  
  wCR2_id <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2")
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", inverse_var = FALSE), wCR2_id)
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n)), wCR2_id)
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = rep(1, n), inverse_var = FALSE), wCR2_id)
  
  wCR2_iv <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
  wCR2_target <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE)
  expect_false(identical(wCR2_target, wCR2_id))
  expect_equal(matrix(wCR2_target, dim(wCR2_target)), matrix(wCR2_iv, dim(wCR2_iv)))
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / dat$w, inverse_var = TRUE), wCR2_target)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = TRUE), CR4_iv)
  
  attr(CR4_iv, "inverse_var") <- FALSE
  CR4_not <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", inverse_var = FALSE)
  expect_equal(CR4_not, CR4_iv)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n)), CR4_not)
  expect_equal(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(lm_fit, cluster = dat$cluster, type = "CR4", target = 1 / dat$w), CR4_not))
  
  wCR4_id <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", inverse_var = FALSE), wCR4_id)
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n)), wCR4_id)
  expect_equal(vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = rep(1, n), inverse_var = FALSE), wCR4_id)
  
  wCR4_iv <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", inverse_var = TRUE)
  wCR4_target <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR4", target = 1 / dat$w, inverse_var = TRUE)
  expect_false(identical(wCR4_target, wCR4_id))
  expect_equal(matrix(wCR4_target, dim(wCR4_target)), matrix(wCR4_iv, dim(wCR4_iv)))
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
  # add tests for CR3f
  CR3f <- vcovCR(lm_fit, cluster = dat$row, type = "CR3f")
  expect_equal(vcovHC(lm_fit, type = "HC3"), nobs(lm_fit) / (nobs(lm_fit) - 1) *as.matrix(CR3f))
  CR3f <- vcovCR(WLS_fit, cluster = dat$row, type = "CR3f")
  expect_equal(vcovHC(WLS_fit, type = "HC3"), nobs(WLS_fit) / (nobs(WLS_fit) - 1) *as.matrix(CR3f))
  
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
  
  df <- m^2 * (m0 - 1) * (m1 - 1) / (m0^2 * (m0 - 1) + m1^2 * (m1 - 1))
  expect_equal(t_Satt$df, df)
})

test_that("Order doesn't matter.",{
  
  check_sort_order(WLS_fit, dat, "cluster")
  
})

test_that("clubSandwich works with dropped observations", {
  dat_miss <- dat
  miss_indicator <- sample.int(n, size = round(n / 10))
  dat_miss$X1[miss_indicator] <- NA
  dat_miss$cluster[miss_indicator] <- NA
  
  lm_dropped <- lm(y ~ X1 + X2 + X3, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(X1))
  lm_complete <- lm(y ~ X1 + X2 + X3, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, cluster = dat_miss$cluster, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, cluster = dat_complete$cluster, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, cluster = dat_miss$cluster, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, cluster = dat_complete$cluster, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


test_that("clubSandwich requires no missing values on the clustering variable", {
  dat_miss <- dat
  miss_indicator <- sample.int(n, size = round(n / 10))
  dat_miss$cluster[miss_indicator] <- NA
  
  lm_dropped <- lm(y ~ X1 + X2 + X3, data = dat_miss)
  
  expect_error(vcovCR(lm_dropped, cluster = dat_miss$cluster, type = "CR0"), 
               "Clustering variable cannot have missing values.")
  expect_error(coef_test(lm_dropped, vcov = "CR0", cluster = dat_miss$cluster, test = "All"),
               "Clustering variable cannot have missing values.")
})


test_that("clubSandwich works with aliased predictors", {
  data(npk, package = "datasets")
  npk_alias <- lm(yield ~ block + N*P*K, npk)
  npk_drop <- lm(yield ~ block + N + P + K + N:P + N:K + P:K, npk)
  
  CR_alias <- lapply(CR_types[-4], function(x) vcovCR(npk_alias, cluster = npk$block, type = x))
  CR_drop <- lapply(CR_types[-4], function(x) vcovCR(npk_drop, cluster = npk$block, type = x))
  expect_equal(CR_alias, CR_drop)
  
  test_drop <- lapply(CR_types[-4], function(x) coef_test(npk_alias, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite"), coefs = 7:12, p_values = FALSE))
  test_complete <- lapply(CR_types[-4], function(x) coef_test(npk_drop, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite"), coefs = 7:12, p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


test_that("weight scale doesn't matter", {
  
  lm_fit_w <- lm(y ~ X1 + X2 + X3, data = dat, weights = rep(4, nrow(dat)))
  
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit, cluster = cluster, type = x))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit_w, cluster = cluster, type = x))
  
  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix), 
               tol = 5 * 10^-7)  
  
  target <- 1 + rpois(nrow(dat), lambda = 8)
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit, cluster = cluster, type = x, target = target))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit_w, cluster = cluster, type = x, target = target * 15))

  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix), 
               tol = 5 * 10^-7)  
  
})

test_that("clubSandwich works with weights of zero.", {
  
  data("LifeCycleSavings")
  n_life <- nrow(LifeCycleSavings)
  LifeCycleSavings$cl <- substr(rownames(LifeCycleSavings), 1, 1)
  table(LifeCycleSavings$cl)
  LifeCycleSavings$wt <- rpois(n_life, lambda = 0.8)
  table(LifeCycleSavings$wt)
  
  lm_full <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings, weights = wt)
  LCS_sub <- subset(LifeCycleSavings, wt > 0)
  lm_sub <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LCS_sub, weights = wt)
  
  CR_full <- lapply(CR_types, function(x) vcovCR(lm_full, cluster = LifeCycleSavings$cl, type = x))
  CR_sub <- lapply(CR_types, function(x) vcovCR(lm_sub, cluster = LCS_sub$cl, type = x))
  expect_equal(CR_full, CR_sub, check.attributes = FALSE)
  
  test_full <- lapply(CR_types, function(x) coef_test(lm_full, vcov = x, cluster = LifeCycleSavings$cl, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  test_sub <- lapply(CR_types, function(x) coef_test(lm_sub, vcov = x, cluster = LCS_sub$cl, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  expect_equal(test_full, test_sub, check.attributes = FALSE)
  
  dat_miss <- LifeCycleSavings
  miss_indicator <- sample.int(n_life, size = round(n_life / 5))
  dat_miss$pop15[miss_indicator] <- NA
  dat_miss$cl[miss_indicator] <- NA
  with(dat_miss, table(wt, is.na(pop15)))
  
  lm_dropped <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = dat_miss, weights = wt)
  dat_complete <- subset(dat_miss, !is.na(pop15))
  lm_complete <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = dat_complete, weights = wt)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, cluster = dat_miss$cl, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, cluster = dat_complete$cl, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, cluster = dat_miss$cl, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, cluster = dat_complete$cl, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
  
})

test_that("vcovCR errors when there is only one cluster.", {
  
  single_cluster_error_msg <- "Cluster-robust variance estimation will not work when the data only includes a single cluster."
  
  expect_error(
    vcovCR(lm_fit, cluster = dat$const, type = "CR0"), 
    single_cluster_error_msg
  )
  
  expect_error(
    conf_int(WLS_fit, vcov = "CR1", cluster = dat$const), single_cluster_error_msg
  )
  
  expect_error(
    coef_test(lm_fit, vcov = "CR2", cluster = dat$const), single_cluster_error_msg
  )
  
  expect_error(
    Wald_test(WLS_fit, constraints = constrain_zero(2:4), 
              vcov = "CR3", cluster = dat$const),
    single_cluster_error_msg
  )
  
})


test_that("vcovCR works with intercept-only model and user-specified weights.", {
  
  lm_int <- lm(y ~ 1, data = dat)
  HC_un <- coef_test(lm_int, vcov="CR2", cluster=dat$row, test = "All")
  
  # Unweighted, HC-robust
  N <- nobs(lm_int)
  yi <- dat$y
  wi <- rep(1, N)
  W <- sum(wi)
  wi <- wi / W
  vi <- rep(1, N)
  V <- sum(vi)
  ei <- residuals_CS(lm_int)
  M <- sum(wi^2 * vi)
  ai <- 1 / sqrt(1 - 2 * wi + M / vi)
  V_hand <- sum(wi^2 * ai^2 * ei^2)
  pi_theta_pj <- diag(vi) - tcrossprod(rep(1,N), wi * vi) - tcrossprod(wi * vi, rep(1, N)) + M
  df <- M^2 / sum(tcrossprod(ai^2 * wi^2) * (pi_theta_pj^2))
  
  expect_true(check_bread(lm_int, cluster = dat$row, y = dat$y))
  expect_true(check_CR(lm_int, vcov = "CR2", cluster = dat$row))
  expect_equal(sqrt(V_hand), HC_un$SE)
  expect_equal(df, HC_un$df_Satt)
  expect_equal(Inf, HC_un$df_z)
  expect_equal(N - 1, HC_un$df_t)
  

  # Unweighted, cluster-robust
  CR_un <- coef_test(lm_int, vcov="CR2", cluster=dat$cluster, test = "All")

  J <- nlevels(dat$cluster)
  w_j <- tapply(wi, dat$cluster, sum)
  e_j <- tapply(wi * ei, dat$cluster, sum) / w_j
  v_j <- tapply(wi^2 * vi, dat$cluster, sum) / w_j^2
  a_j <- 1 / sqrt(1 - 2 * w_j + M / v_j)
  V_hand <- sum(w_j^2 * a_j^2 * e_j^2)
  pi_theta_pj <- diag(v_j) - tcrossprod(rep(1,J), w_j * v_j) - tcrossprod(w_j * v_j, rep(1, J)) + M
  df <- M^2 / sum(tcrossprod(a_j^2 * w_j^2) * (pi_theta_pj^2))
  
  expect_true(check_bread(lm_int, cluster = dat$cluster, y = dat$y))
  expect_true(check_CR(lm_int, vcov = "CR2", cluster = dat$cluster))
  expect_equal(sqrt(V_hand), CR_un$SE)
  expect_equal(df, CR_un$df_Satt)
  expect_equal(Inf, CR_un$df_z)
  expect_equal(J - 1, CR_un$df_t)

  
  # Weighted, HC-robust  
  WLS_int <- lm(y ~ 1, data = dat, weights = w)
  HC_wt <- coef_test(WLS_int, vcov="CR2", cluster=dat$row, test = "All")
  
  N <- nobs(WLS_int)
  yi <- dat$y
  wi <- WLS_int$weights
  W <- sum(wi)
  wi <- wi / W
  vi <- rep(1, N)
  V <- sum(vi)
  ei <- residuals_CS(WLS_int)
  M <- sum(wi^2 * vi)
  ai <- 1 / sqrt(1 - 2 * wi + M / vi)
  V_hand <- sum(wi^2 * ai^2 * ei^2)
  pi_theta_pj <- diag(vi) - tcrossprod(rep(1,N), wi * vi) - tcrossprod(wi * vi, rep(1, N)) + M
  df <- M^2 / sum(tcrossprod(ai^2 * wi^2) * (pi_theta_pj^2))
  
  expect_true(check_bread(WLS_int, cluster = dat$row, y = dat$y))
  expect_true(check_CR(WLS_int, vcov = "CR2", cluster = dat$row))
  expect_equal(sqrt(V_hand), HC_wt$SE)
  expect_equal(df, HC_wt$df_Satt)
  expect_equal(Inf, HC_wt$df_z)
  expect_equal(N - 1, HC_wt$df_t)
  
  # Weighted, cluster-robust
  CR_wt <- coef_test(WLS_int, vcov="CR2", cluster=dat$cluster, test = "All")
  
  J <- nlevels(dat$cluster)
  w_j <- tapply(wi, dat$cluster, sum)
  e_j <- tapply(wi * ei, dat$cluster, sum) / w_j
  v_j <- tapply(wi^2 * vi, dat$cluster, sum) / w_j^2
  a_j <- 1 / sqrt(1 - 2 * w_j + M / v_j)
  V_hand <- sum(w_j^2 * a_j^2 * e_j^2)
  pi_theta_pj <- diag(v_j) - tcrossprod(rep(1,J), w_j * v_j) - tcrossprod(w_j * v_j, rep(1, J)) + M
  df <- M^2 / sum(tcrossprod(a_j^2 * w_j^2) * (pi_theta_pj^2))
  
  expect_true(check_bread(WLS_int, cluster = dat$cluster, y = dat$y))
  expect_true(check_CR(WLS_int, vcov = "CR2", cluster = dat$cluster))
  expect_equal(sqrt(V_hand), CR_wt$SE, tolerance = 10^-3)
  expect_equal(df, CR_wt$df_Satt, tolerance = 10^-3)
  expect_equal(Inf, CR_wt$df_z)
  expect_equal(J - 1, CR_wt$df_t)
  
})
