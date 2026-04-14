context("estimatr::iv_robust objects")

skip_if_not_installed("estimatr")

suppressMessages(library(estimatr, quietly = TRUE))

set.seed(20190513)

# -------------------------------------------------
# Set up test data with IV structure
# -------------------------------------------------

N <- 200
J <- 20
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

# Build manual matrices for verification
mf <- model.frame(iv_un$terms, data = dat)
X <- model.matrix(iv_un$terms_regressors, data = mf)
inst_formula <- as.formula(paste("~", deparse(iv_un$formula[[3]][[3]])))
Z <- model.matrix(inst_formula, data = mf)
y_vec <- mf$y


# -------------------------------------------------
# Tests: model.frame
# -------------------------------------------------

test_that("model.frame() works for iv_robust", {

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
})


# -------------------------------------------------
# Tests: model_matrix (projected X)
# -------------------------------------------------

test_that("model_matrix() returns projected matrix with correct dimensions", {

  mm <- model_matrix(iv_un)
  expect_equal(nrow(mm), N)
  expect_equal(ncol(mm), ncol(X))  # same number of columns as regressors
  expect_equal(colnames(mm), colnames(X))
})

test_that("model_matrix() agrees with manual projection for unweighted model", {

  XZ <- model_matrix(iv_un)
  ZtZ_inv <- chol2inv(chol(crossprod(Z)))
  XZ_check <- Z %*% ZtZ_inv %*% crossprod(Z, X)

  expect_equal(XZ, XZ_check, check.attributes = FALSE)
})

test_that("model_matrix() agrees with manual projection for weighted model", {

  XZ_wt <- model_matrix(iv_wt)
  ZwZ_inv <- chol2inv(chol(crossprod(Z, w * Z)))
  XZ_check <- Z %*% ZwZ_inv %*% crossprod(Z, w * X)

  expect_equal(XZ_wt, XZ_check, check.attributes = FALSE)
})


# -------------------------------------------------
# Tests: residuals
# -------------------------------------------------

test_that("residuals_CS() matches y - X * beta", {

  r_un <- residuals_CS(iv_un)
  r_check <- as.vector(y_vec - X %*% coef(iv_un))
  expect_equal(as.vector(r_un), r_check)

  r_wt <- residuals_CS(iv_wt)
  r_check_wt <- as.vector(y_vec - X %*% coef(iv_wt))
  expect_equal(as.vector(r_wt), r_check_wt)
})


# -------------------------------------------------
# Tests: bread
# -------------------------------------------------

test_that("bread() agrees with manual computation for unweighted model", {

  XZ <- model_matrix(iv_un)
  bread_check <- chol2inv(chol(crossprod(XZ))) * nobs(iv_un)
  expect_equal(bread(iv_un), bread_check, check.attributes = FALSE)
})

test_that("bread() agrees with manual computation for weighted model", {

  XZ <- model_matrix(iv_wt)
  bread_check <- chol2inv(chol(crossprod(XZ, w * XZ))) * nobs(iv_wt)
  expect_equal(bread(iv_wt), bread_check, check.attributes = FALSE)
})

test_that("bread works with check_bread()", {

  expect_true(check_bread(iv_un, cluster = dat$cluster, y = dat$y))
  expect_true(check_bread(iv_wt, cluster = dat$cluster, y = dat$y))
})


# -------------------------------------------------
# Tests: vcovCR runs with explicit args
# -------------------------------------------------

test_that("vcovCR works with explicit cluster and type", {

  for (cr in CR_types) {
    V <- vcovCR(iv_un, cluster = dat$cluster, type = cr)
    expect_s3_class(V, "clubSandwich")
    expect_equal(nrow(V), length(coef(iv_un)))
    expect_equal(ncol(V), length(coef(iv_un)))
  }
})

test_that("vcovCR works for weighted model", {

  for (cr in CR_types) {
    V <- vcovCR(iv_wt, cluster = dat$cluster, type = cr)
    expect_s3_class(V, "clubSandwich")
  }
})


# -------------------------------------------------
# Tests: cluster auto-detection
# -------------------------------------------------

test_that("vcovCR auto-detects cluster from iv_robust object", {

  V_auto <- vcovCR(iv_cl)
  V_explicit <- vcovCR(iv_cl, cluster = dat$cluster, type = "CR2")
  expect_equal(V_auto, V_explicit, check.attributes = FALSE)
})

test_that("vcovCR errors when cluster is missing on unclustered model", {

  expect_error(
    vcovCR(iv_un, type = "CR2"),
    "You must specify a clustering variable"
  )
})


# -------------------------------------------------
# Tests: SE type inference
# -------------------------------------------------

test_that("vcovCR inherits se_type from iv_robust object", {

  # CR2
  V_auto_cr2 <- vcovCR(iv_cl)
  V_explicit_cr2 <- vcovCR(iv_cl, cluster = dat$cluster, type = "CR2")
  expect_equal(V_auto_cr2, V_explicit_cr2, check.attributes = FALSE)

  # CR0
  V_auto_cr0 <- vcovCR(iv_cl_cr0)
  V_explicit_cr0 <- vcovCR(iv_cl_cr0, cluster = dat$cluster, type = "CR0")
  expect_equal(V_auto_cr0, V_explicit_cr0, check.attributes = FALSE)

  # stata -> CR1S
  V_auto_stata <- vcovCR(iv_cl_stata)
  V_explicit_stata <- vcovCR(iv_cl_stata, cluster = dat$cluster, type = "CR1S")
  expect_equal(V_auto_stata, V_explicit_stata, check.attributes = FALSE)
})

test_that("vcovCR errors when se_type cannot be mapped", {

  expect_error(
    vcovCR(iv_un, cluster = dat$cluster),
    "You must specify a `type`"
  )
})


# -------------------------------------------------
# Tests: inverse_var rejection
# -------------------------------------------------

test_that("vcovCR rejects inverse_var = TRUE", {

  expect_error(
    vcovCR(iv_un, cluster = dat$cluster, type = "CR2", inverse_var = TRUE),
    "inverse_var option is not available"
  )
})


# -------------------------------------------------
# Tests: na.action
# -------------------------------------------------

test_that("na.action() works correctly", {

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

test_that("vcovCR options don't matter for CR0", {

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

test_that("vcovCR options work for CR2", {

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

test_that("vcovCR options work for CR4", {

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

test_that("CR2 is target-unbiased", {

  expect_true(check_CR(iv_un, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(iv_wt, vcov = "CR2", cluster = dat$cluster))
})
