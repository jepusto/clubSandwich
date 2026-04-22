context("estimatr::iv_robust objects")

skip_if_not_installed("estimatr")

suppressMessages(library(estimatr, quietly = TRUE))

set.seed(20190513)

# -------------------------------------------------
# Set up test data with IV structure
# -------------------------------------------------

N <- 200L
J <- 20L
cluster <- rep(1:J, each = N / J)

# Instruments
z1 <- rnorm(N)
z2 <- rnorm(N)

# Exogenous regressor
x_exog <- rnorm(N)

# Endogenous regressor (correlated with instruments)
x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(N)

# Outcome
y <- 1 + 2 * x_endog + 0.5 * x_exog + rnorm(N)

# Weights
w <- 1 + rpois(N, lambda = 3)

dat <- data.frame(
  y = y, x_endog = x_endog, x_exog = x_exog,
  z1 = z1, z2 = z2, w = w, cluster = cluster
)

CR_types <- paste0("CR", 0:4)

# -------------------------------------------------
# Fit models
# -------------------------------------------------

# Unweighted, no clusters in call
iv_un <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat
)

# Weighted, no clusters in call
iv_wt <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, weights = w
)

# Clustered with CR2
iv_cl <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, clusters = cluster, se_type = "CR2"
)

# Clustered with CR0
iv_cl_cr0 <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, clusters = cluster, se_type = "CR0"
)

# Clustered with stata
iv_cl_stata <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, clusters = cluster, se_type = "stata"
)

# Weighted, no clusters in call
iv_cl_wt <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, 
  weights = w,
  clusters = cluster
)


# Build manual matrices for verification
mf <- model.frame(iv_un$terms, data = dat)
X <- model.matrix(iv_un$terms_regressors, data = mf)
inst_formula <- as.formula(paste("~", deparse(iv_un$formula[[3]][[3]])))
Z <- model.matrix(inst_formula, data = mf)
y_vec <- mf$y


# -------------------------------------------------
# Tests: model.frame
# -------------------------------------------------

test_that("model.frame() works for iv_robust.", {

  mf_un <- model.frame(iv_un)
  expect_equal(ncol(mf_un), 5L)  # y, x_endog, x_exog, z1, z2
  expect_true(all(c("y", "x_endog", "x_exog", "z1", "z2") %in% names(mf_un)))
  expect_equal(nrow(mf_un), N)

  # Clustered model has (clusters) column
  mf_cl <- model.frame(iv_cl)
  expect_true("(clusters)" %in% names(mf_cl))
  expect_equal(mf_cl[["(clusters)"]], dat$cluster)

  # Weighted model has (weights) column
  mf_wt <- model.frame(iv_wt)
  expect_true("(weights)" %in% names(mf_wt))

  # Clustered and weighted model has (clusters) and (weights) columns
  mf_cl_wt <- model.frame(iv_cl_wt)
  expect_true("(weights)" %in% names(mf_cl_wt))
  expect_true("(clusters)" %in% names(mf_cl_wt))
  expect_equal(mf_cl[["(clusters)"]], dat$cluster)
  
})


# -------------------------------------------------
# Tests: model_matrix (projected X)
# -------------------------------------------------

test_that("model_matrix() returns projected matrix with correct dimensions.", {

  mm <- model_matrix(iv_un)
  expect_identical(dim(mm), c(N,ncol(X)))
  expect_equal(colnames(mm), colnames(X))
})

test_that("model_matrix() agrees with manual projection.", {

  # Unweighted model
  XZ <- model_matrix(iv_un)
  ZtZ_inv <- chol2inv(chol(crossprod(Z)))
  XZ_check <- Z %*% ZtZ_inv %*% crossprod(Z, X)

  expect_equal(XZ, XZ_check, check.attributes = FALSE)
  
  # Weighted model
  XZ_wt <- model_matrix(iv_wt)
  ZwZ_inv <- chol2inv(chol(crossprod(Z, w * Z)))
  XZ_check <- Z %*% ZwZ_inv %*% crossprod(Z, w * X)
  
  expect_equal(XZ_wt, XZ_check, check.attributes = FALSE)
  
})

# -------------------------------------------------
# Tests: residuals
# -------------------------------------------------

test_that("residuals_CS() matches y - X * beta.", {

  r_un <- residuals_CS(iv_un)
  r_check <- as.vector(y_vec - X %*% coef_CS(iv_un))
  expect_equal(as.vector(r_un), r_check)

  r_wt <- residuals_CS(iv_wt)
  r_check_wt <- as.vector(y_vec - X %*% coef_CS(iv_wt))
  expect_equal(as.vector(r_wt), r_check_wt)
})


# -------------------------------------------------
# Tests: bread
# -------------------------------------------------

test_that("bread() agrees with manual computation.", {

  # Unweighted model
  XZ <- model_matrix(iv_un)
  bread_check <- chol2inv(chol(crossprod(XZ))) * nobs(iv_un)
  expect_equal(bread(iv_un), bread_check, check.attributes = FALSE)
  expect_true(check_bread(iv_un, cluster = dat$cluster, y = dat$y))
  
  # Weighted model
  XZ <- model_matrix(iv_wt)
  bread_check <- chol2inv(chol(crossprod(XZ, w * XZ))) * nobs(iv_wt)
  expect_equal(bread(iv_wt), bread_check, check.attributes = FALSE)
  expect_true(check_bread(iv_wt, cluster = dat$cluster, y = dat$y))
})

# -------------------------------------------------
# Tests: vcovCR runs with explicit args
# -------------------------------------------------

test_that("vcovCR works with explicit cluster and type.", {

  for (cr in CR_types) {
    V <- vcovCR(iv_un, cluster = dat$cluster, type = cr)
    expect_s3_class(V, "clubSandwich")
    expect_identical(dim(V), rep(length(coef(iv_un)), 2L))
  }
})

test_that("vcovCR works for weighted model.", {

  for (cr in CR_types) {
    V <- vcovCR(iv_wt, cluster = dat$cluster, type = cr)
    expect_s3_class(V, "clubSandwich")
  }
})


# -------------------------------------------------
# Tests: cluster auto-detection
# -------------------------------------------------

test_that("vcovCR auto-detects cluster from iv_robust object.", {

  V_auto <- vcovCR(iv_cl)
  V_explicit <- vcovCR(iv_cl, cluster = dat$cluster, type = "CR2")
  expect_equal(V_auto, V_explicit, check.attributes = FALSE)
})

test_that("vcovCR errors when cluster is missing on unclustered model.", {

  expect_error(
    vcovCR(iv_un, type = "CR2"),
    "You must specify a clustering variable"
  )
})


# -------------------------------------------------
# Tests: SE type inference
# -------------------------------------------------

test_that("vcovCR inherits se_type from iv_robust object.", {

  # CR2
  V_auto_cr2 <- vcovCR(iv_cl)
  V_explicit_cr2 <- vcovCR(iv_cl, cluster = dat$cluster, type = "CR2")
  expect_equal(V_auto_cr2, V_explicit_cr2, check.attributes = FALSE)

  V_auto_wt_cr2 <- vcovCR(iv_cl_wt)
  V_explicit_wt_cr2 <- vcovCR(iv_cl_wt, cluster = dat$cluster, type = "CR2")
  expect_equal(V_auto_wt_cr2, V_explicit_wt_cr2, check.attributes = FALSE)
  
  # CR0
  V_auto_cr0 <- vcovCR(iv_cl_cr0)
  V_explicit_cr0 <- vcovCR(iv_cl_cr0, cluster = dat$cluster, type = "CR0")
  expect_equal(V_auto_cr0, V_explicit_cr0, check.attributes = FALSE)

  # stata -> CR1S
  V_auto_stata <- vcovCR(iv_cl_stata)
  V_explicit_stata <- vcovCR(iv_cl_stata, cluster = dat$cluster, type = "CR1S")
  expect_equal(V_auto_stata, V_explicit_stata, check.attributes = FALSE)
})

test_that("vcovCR errors when se_type cannot be mapped.", {

  expect_error(
    vcovCR(iv_un, cluster = dat$cluster),
    "You must specify a `type`"
  )
})


# -------------------------------------------------
# Tests: inverse_var rejection
# -------------------------------------------------

test_that("vcovCR rejects inverse_var = TRUE.", {

  expect_error(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR2", inverse_var = TRUE),
    "inverse_var option is not available"
  )
})


# -------------------------------------------------
# Tests: na.action
# -------------------------------------------------

test_that("na.action() works correctly.", {

  dat_na <- dat
  dat_na$y[c(5, 15, 25)] <- NA
  iv_na <- iv_robust(
    y ~ x_endog + x_exog | x_exog + z1 + z2,
    data = dat_na, clusters = cluster, se_type = "CR2"
  )

  na_act <- na.action(iv_na)
  expect_s3_class(na_act, "omit")
  expect_equal(as.integer(na_act), c(5L, 15L, 25L))
  expect_equal(nobs(iv_na), N - 3L)
})


# -------------------------------------------------
# Tests: vcovCR options
# -------------------------------------------------

test_that("vcovCR options don't matter for CR0.", {

  CR0 <- vcovCR(iv_un, cluster = dat$cluster, type = "CR0")
  expect_output(print(CR0))

  attr(CR0, "target") <- NULL
  attr(CR0, "inverse_var") <- NULL

  CR0_A <- vcovCR(iv_un, cluster = dat$cluster, type = "CR0", target = 1 / dat$w)
  attr(CR0_A, "target") <- NULL
  attr(CR0_A, "inverse_var") <- NULL
  expect_identical(CR0_A, CR0)

  CR0_B <- vcovCR(iv_un, cluster = dat$cluster, type = "CR0",
                   target = 1 / dat$w, inverse_var = FALSE)
  attr(CR0_B, "target") <- NULL
  attr(CR0_B, "inverse_var") <- NULL
  expect_identical(CR0_B, CR0)

  expect_error(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR0",
           target = 1 / dat$w, inverse_var = TRUE)
  )
})

test_that("vcovCR options work for CR2.", {

  CR2 <- vcovCR(iv_un, cluster = dat$cluster, type = "CR2")
  expect_equal(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR2",
           target = rep(1, nobs(iv_un))),
    CR2
  )
  expect_false(identical(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR2",
           target = 1 / dat$w),
    CR2
  ))
})

test_that("vcovCR options work for CR4.", {

  CR4 <- vcovCR(iv_un, cluster = dat$cluster, type = "CR4")
  expect_identical(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR4",
           target = rep(1, nobs(iv_un))),
    CR4
  )
  expect_false(identical(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR4",
           target = 1 / dat$w),
    CR4
  ))
})


# -------------------------------------------------
# Tests: CR2 target-unbiasedness
# -------------------------------------------------

test_that("CR2 is target-unbiased.", {

  expect_true(check_CR(iv_un, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(iv_wt, vcov = "CR2", cluster = dat$cluster))
})

test_that("CR4 is target-unbiased.", {
  skip("Need to understand target-unbiasedness for iv_robust objects.")
  expect_true(check_CR(iv_un, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(iv_wt, vcov = "CR4", cluster = dat$cluster))
})


# -------------------------------------------------
# Tests: equivalence to vcovHC with size-1 clusters
# -------------------------------------------------

test_that("vcovCR is equivalent to vcovHC (with HC0-HC3) when clusters are all of size 1.", {
  suppressMessages(library(sandwich, quietly = TRUE))
  iv_fit <- iv_robust(
    y ~ x_endog + x_exog | x_exog + z1 + z2,
    data = dat, se_type = "HC0"
  )
  CR0 <- vcovCR(iv_fit, cluster = 1:nobs(iv_fit), type = "CR0")
  expect_equal(as.matrix(CR0), vcov(iv_fit), check.attributes = FALSE)
  
  CR1 <- vcovCR(iv_fit, cluster = 1:nobs(iv_fit), type = "CR1p")
  expect_equal(as.matrix(CR1), vcov(update(iv_fit, se_type = "HC1")), check.attributes = FALSE)
  
  CR2 <- vcovCR(iv_fit, cluster = 1:nobs(iv_fit), type = "CR2")
  expect_equal(as.matrix(CR2), vcov(update(iv_fit, se_type = "HC2")), check.attributes = FALSE)

  CR3 <- vcovCR(iv_fit, cluster = 1:nobs(iv_fit), type = "CR3")
  expect_equal(as.matrix(CR3), vcov(update(iv_fit, se_type = "HC3")), check.attributes = FALSE)
  
})


# -------------------------------------------------
# Tests: sort-order invariance
# -------------------------------------------------

test_that("Order doesn't matter.", {

  check_sort_order(iv_un, dat, "cluster")
  check_sort_order(iv_wt, dat, "cluster")
  check_sort_order(iv_cl, dat, "cluster")
  check_sort_order(iv_cl_wt, dat, "cluster")
  
})


# -------------------------------------------------
# Tests: dropped observations
# -------------------------------------------------

test_that("clubSandwich works with dropped observations", {
  dat_miss <- dat
  dat_miss$x_exog[sample.int(N, size = round(N / 10))] <- NA
  iv_dropped <- update(iv_un, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(x_exog))
  iv_complete <- update(iv_un, data = dat_complete)

  CR_drop <- lapply(CR_types, function(x) vcovCR(iv_dropped, cluster = dat_miss$cluster, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(iv_complete, cluster = dat_complete$cluster, type = x))
  expect_equal(CR_drop, CR_complete)

  test_drop <- lapply(CR_types, function(x) coef_test(iv_dropped, vcov = x, cluster = dat_miss$cluster, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(iv_complete, vcov = x, cluster = dat_complete$cluster, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


# -------------------------------------------------
# Tests: weight scale invariance
# -------------------------------------------------

test_that("weight scale doesn't matter", {

  iv_fit_w <- update(iv_un, weights = rep(4, nobs(iv_un)))

  unweighted_fit <- lapply(CR_types, function(x) vcovCR(iv_un, cluster = dat$cluster, type = x))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(iv_fit_w, cluster = dat$cluster, type = x))

  expect_equal(lapply(unweighted_fit, as.matrix),
               lapply(weighted_fit, as.matrix),
               tol = 5 * 10^-7)

  target <- 1 + rpois(N, lambda = 8)
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(iv_un, cluster = dat$cluster, type = x, target = target))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(iv_fit_w, cluster = dat$cluster, type = x, target = target * 15))

  expect_equal(lapply(unweighted_fit, as.matrix),
               lapply(weighted_fit, as.matrix),
               tol = 5 * 10^-7)

})


# -------------------------------------------------
# Tests: weights of zero
# -------------------------------------------------

test_that("clubSandwich works with weights of zero.", {

  dat$awt <- rpois(N, lambda = 1.4)

  iv_full <- update(iv_un, weights = awt)
  dat_sub <- subset(dat, awt > 0)
  iv_sub <- update(iv_full, data = dat_sub)

  CR_full <- lapply(CR_types, function(x) vcovCR(iv_full, cluster = dat$cluster, type = x))
  CR_sub <- lapply(CR_types, function(x) vcovCR(iv_sub, cluster = dat_sub$cluster, type = x))
  expect_equal(CR_full, CR_sub, check.attributes = FALSE)

  test_full <- lapply(CR_types, function(x) coef_test(iv_full, vcov = x, cluster = dat$cluster, test = c("z", "naive-t", "Satterthwaite"), p_values = TRUE))
  test_sub <- lapply(CR_types, function(x) coef_test(iv_sub, vcov = x, cluster = dat_sub$cluster, test = c("z", "naive-t", "Satterthwaite"), p_values = TRUE))
  expect_equal(test_full, test_sub, check.attributes = FALSE)

})


# -------------------------------------------------
# Tests: fixed_effects
# -------------------------------------------------

iv_fe <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat, clusters = cluster,
  fixed_effects = ~ cluster, se_type = "CR2"
)

test_that("model.frame() includes fixed-effects variables", {
  mf_fe <- model.frame(iv_fe)
  expect_true("cluster" %in% names(mf_fe))
})

test_that("augmented_model_matrix() returns the FE design matrix", {
  amm <- augmented_model_matrix(iv_fe, cluster = dat$cluster,
                                inverse_var = FALSE, ignore_FE = FALSE)
  expect_equal(nrow(amm), N)
  expect_equal(ncol(amm), J)
  expect_equal(colnames(amm), paste0("cluster", levels(as.factor(dat$cluster))))

  # For a non-FE model, augmented_model_matrix() returns NULL
  expect_null(augmented_model_matrix(iv_un, cluster = dat$cluster,
                                     inverse_var = FALSE, ignore_FE = FALSE))
})

test_that("model_matrix() returns FE-residualized projected matrix", {
  XP <- model_matrix(iv_fe)
  expect_equal(nrow(XP), N)
  expect_equal(ncol(XP), length(coef(iv_fe)))

  # Coefficients should be recoverable from X_proj and FE-residualized outcome
  F_mat <- model.matrix(~ 0 + cluster, data = dat)
  y_demean <- lm.fit(F_mat, dat$y)$residuals
  expect_equal(lm.fit(XP, y_demean)$coefficients, coef(iv_fe),
               check.attributes = FALSE)

  # Projected columns should be centered within each cluster
  within_means <- apply(XP, 2, function(x) tapply(x, dat$cluster, mean))
  expect_lt(max(abs(within_means)), 1e-12)
})

test_that("model_matrix() FE path agrees without fixest", {
  local_mocked_bindings(
    requireNamespace = function(...) FALSE
  )
  expect_false(requireNamespace("fixest", quietly = TRUE))

  XP <- model_matrix(iv_fe)
  F_mat <- model.matrix(~ 0 + cluster, data = dat)
  y_demean <- lm.fit(F_mat, dat$y)$residuals
  expect_equal(lm.fit(XP, y_demean)$coefficients, coef(iv_fe),
               check.attributes = FALSE)

  expect_message(
    vcovCR(iv_fe, type = "CR0"),
    "For improved performance in models with fixed effects, install the \\{fixest\\} package."
  )
})

test_that("vcovCR works with fixed effects", {

  for (cr in CR_types) {
    V <- vcovCR(iv_fe, cluster = dat$cluster, type = cr)
    expect_s3_class(V, "clubSandwich")
    expect_equal(nrow(V), length(coef(iv_fe)))
    expect_equal(ncol(V), length(coef(iv_fe)))
  }
})

test_that("bread() works with fixed effects", {
  expect_true(check_bread(iv_fe, cluster = dat$cluster, y = dat$y))
})

test_that("CR2 is target-unbiased with fixed effects", {
  expect_true(check_CR(iv_fe, vcov = "CR2"))
})

test_that("vcovCR matches iv_robust$vcov for FE model with se_type='CR0'", {
  iv_fe_cr0 <- iv_robust(
    y ~ x_endog + x_exog | x_exog + z1 + z2,
    data = dat, clusters = cluster,
    fixed_effects = ~ cluster, se_type = "CR0"
  )
  V_cr0 <- vcovCR(iv_fe_cr0, type = "CR0")
  expect_equal(as.matrix(V_cr0), vcov(iv_fe_cr0), check.attributes = FALSE)
})


# -------------------------------------------------
# Tests: missing cluster
# -------------------------------------------------

test_that("clubSandwich requires no missing values on the clustering variable", {

  dat_miss <- dat
  miss_indicator <- sample.int(N, size = round(N / 10))
  dat_miss$cluster[miss_indicator] <- NA

  iv_dropped <- iv_robust(y ~ x_endog + x_exog | x_exog + z1 + z2, data = dat_miss)
  expect_error(vcovCR(iv_dropped, cluster = dat_miss$cluster, type = "CR0"),
               "Clustering variable cannot have missing values.")
  expect_error(coef_test(iv_dropped, vcov = "CR0", cluster = dat_miss$cluster, test = "All"),
               "Clustering variable cannot have missing values.")
})


# -------------------------------------------------
# Higher-level inference tests: compare against ivreg::ivreg
# -------------------------------------------------

skip_if_not_installed("ivreg")

iv_reg <- ivreg::ivreg(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat
)
iv_est <- iv_robust(
  y ~ x_endog + x_exog | x_exog + z1 + z2,
  data = dat
)

test_that("Wald_test() agrees between iv_robust and ivreg", {

  W_reg <- Wald_test(iv_reg, constraints = constrain_zero("x_endog"),
                     vcov = "CR2", cluster = dat$cluster)
  W_est <- Wald_test(iv_est, constraints = constrain_zero("x_endog"),
                     vcov = "CR2", cluster = dat$cluster)
  expect_equal(W_reg, W_est)

})


test_that("conf_int() agrees between iv_robust and ivreg", {

  ci_reg <- conf_int(iv_reg, vcov = "CR2", cluster = dat$cluster)
  ci_est <- conf_int(iv_est, vcov = "CR2", cluster = dat$cluster)
  expect_equal(ci_reg, ci_est)

})


test_that("coef_test() agrees between iv_robust and ivreg", {

  ct_reg <- coef_test(iv_reg, vcov = "CR2", cluster = dat$cluster)
  ct_est <- coef_test(iv_est, vcov = "CR2", cluster = dat$cluster)
  expect_equal(ct_reg, ct_est)

})
