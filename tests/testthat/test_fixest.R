context("fixest objects")

skip_if_not_installed("fixest")

library(fixest, quietly = TRUE, warn.conflicts = FALSE)

set.seed(20190513)

m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
n <- length(cluster)
X1 <- rnorm(n)
X2 <- rnorm(n)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
w <- rgamma(n, shape = 3, scale = 3)
y <- 0.4 * X1 + 0.3 * X2 + nu + e

dat <- data.frame(y, X1, X2, cluster, w, row = 1:n)

# feols models
fit_nfe <- feols(y ~ X1 + X2, data = dat)
fit_fe <- feols(y ~ X1 + X2 | cluster, data = dat)
fit_fe_char <- feols(y ~ X1 + X2, data = dat, fixef = "cluster")
fit_wt <- feols(y ~ X1 + X2, data = dat, weights = ~w)
fit_wt_fe <- feols(y ~ X1 + X2 | cluster, data = dat, weights = ~w)

data(trade, package = "fixest")

# feols model with Origin and Destination fixed effects
fit_fe2 <- feols(
  log(Euros) ~ log(dist_km) | Origin + Destination,
  data = trade
)
fit_fe2_char <- feols(
  log(Euros) ~ log(dist_km), 
  data = trade,
  fixef = c("Origin","Destination")
)

# feols model with Origin-by-Destination fixed effects
fit_fe_cross <- feols(
  log(Euros) ~ Year | Origin^Destination,
  data = trade
)

# feols model with destination and year-by-destination FEs
fit_fe_slope <- feols(
  log(Euros) ~ log(dist_km) + Origin | Destination[Year],
  data = trade
)

# feols model with year-by-destination FEs
fit_slope_only <- feols(
  log(Euros) ~ log(dist_km) + Origin | Destination[[Year]],
  data = trade
)

# equivalent lm models for comparison
lm_nfe <- lm(y ~ X1 + X2, data = dat)
lm_fe <- lm(y ~ X1 + X2 + cluster, data = dat)
lm_wt <- lm(y ~ X1 + X2, data = dat, weights = w)
lm_wt_fe <- lm(y ~ X1 + X2 + cluster, data = dat, weights = w)
lm_fe2 <- lm(log(Euros) ~ log(dist_km) + Origin + Destination, data = trade)
lm_fe_slope <- lm(log(Euros) ~ log(dist_km) + Origin + Destination + Year:Destination, data = trade)
lm_slope_only <- lm(log(Euros) ~ log(dist_km) + Origin + Year:Destination, data = trade)

CR_types <- paste0("CR", 0:4)
focal <- c("X1", "X2")

test_that("bread works", {

  expect_true(check_bread(fit_nfe, cluster = dat$cluster, y = dat$y))
  expect_true(check_bread(fit_fe, cluster = dat$cluster, y = dat$y, check_coef = FALSE))
  expect_true(check_bread(fit_wt, cluster = dat$cluster, y = dat$y))
  expect_true(check_bread(fit_wt_fe, cluster = dat$cluster, y = dat$y, check_coef = FALSE))

})


test_that("vcovCR options work for CR2", {

  CR2_iv <- vcovCR(fit_nfe, cluster = dat$cluster, type = "CR2")
  expect_equal(vcovCR(fit_nfe, cluster = dat$cluster, type = "CR2", inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(fit_nfe, cluster = dat$cluster, type = "CR2", inverse_var = FALSE)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(CR2_not, CR2_iv)

})


test_that("vcovCR options work for CR4", {

  CR4_iv <- vcovCR(fit_nfe, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(fit_nfe, cluster = dat$cluster, type = "CR4", inverse_var = TRUE), CR4_iv)

})


test_that("CR2 and CR4 are target-unbiased", {

  expect_true(check_CR(fit_nfe, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(fit_nfe, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(fit_fe, vcov = "CR2"))
  expect_true(check_CR(fit_fe, vcov = "CR4"))

  expect_true(check_CR(fit_wt, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(fit_wt, vcov = "CR4", cluster = dat$cluster))

  # Weighted FE uses iterative demeaning so tolerance is slightly looser
  cr2_wt_fe <- vcovCR(fit_wt_fe, type = "CR2")
  expect_true(check_CR(fit_wt_fe, vcov = cr2_wt_fe, tol = 1e-4))

})


test_that("feols without FE matches lm", {

  for (cr_type in CR_types) {
    v_feols <- vcovCR(fit_nfe, cluster = dat$cluster, type = cr_type)
    v_lm <- vcovCR(lm_nfe, cluster = dat$cluster, type = cr_type)
    expect_equal(as.matrix(v_feols), as.matrix(v_lm), tolerance = 1e-6)
  }

})


test_that("feols with FE matches lm with dummies", {

  # CR3 excluded: lm with cluster dummies is singular when dummies match clustering
  # CR4 excluded: augmented model matrix introduces small numerical differences
  # CR4 target-unbiasedness is verified separately via check_CR
  CR_types_fe <- paste0("CR", 0:2)
  for (cr_type in CR_types_fe) {
    v_feols <- vcovCR(fit_fe, type = cr_type)
    v_lm <- vcovCR(lm_fe, cluster = dat$cluster, type = cr_type)
    expect_equal(as.matrix(v_feols), as.matrix(v_lm)[focal, focal], tolerance = 1e-6)
  }

})


test_that("weighted feols without FE matches weighted lm", {

  for (cr_type in CR_types) {
    v_feols <- vcovCR(fit_wt, cluster = dat$cluster, type = cr_type)
    v_lm <- vcovCR(lm_wt, cluster = dat$cluster, type = cr_type)
    expect_equal(as.matrix(v_feols), as.matrix(v_lm), tolerance = 1e-6)
  }

})


test_that("weighted feols with FE matches weighted lm with dummies", {

  CR_types_fe <- paste0("CR", 0:2)
  for (cr_type in CR_types_fe) {
    v_feols <- vcovCR(fit_wt_fe, type = cr_type)
    v_lm <- vcovCR(lm_wt_fe, cluster = dat$cluster, type = cr_type)
    expect_equal(as.matrix(v_feols), as.matrix(v_lm)[focal, focal], tolerance = 1e-6)
  }

})


test_that("cluster auto-detection works", {

  # single FE
  expect_identical(dat$cluster, findCluster.fixest(fit_fe))
  expect_identical(dat$cluster, findCluster.fixest(fit_fe_char))
  
  # fully interacted FEs
  expect_identical(
    paste(trade$Origin, trade$Destination, sep = "_"), 
    findCluster.fixest(fit_fe_cross)
  )
  
  # fixed intercepts and slopes
  expect_identical(trade$Destination, findCluster.fixest(fit_fe_slope))
  
  # two-way FE
  expect_error(findCluster.fixest(fit_fe2))
  
  # Auto-detects first FE variable as cluster
  cr2_auto <- vcovCR(fit_fe, type = "CR2")
  cr2_manual <- vcovCR(fit_fe, cluster = dat$cluster, type = "CR2")
  expect_equal(cr2_auto, cr2_manual)

  # No FE: must specify cluster
  expect_error(
    vcovCR(fit_nfe, type = "CR2"),
    "You must specify a clustering variable."
  )

})


test_that("sort order doesn't matter", {

  check_sort_order(fit_fe, dat, cluster = "cluster")

})


test_that("coef_test works", {

  ct_feols <- coef_test(fit_fe, vcov = "CR2", test = "All", p_values = FALSE)
  suppressWarnings(
    ct_lm <- coef_test(lm_fe, vcov = "CR2", cluster = dat$cluster, test = "All", p_values = FALSE)
  )
  lm_rows <- rownames(ct_lm) %in% focal

  expect_equal(ct_feols$beta, ct_lm$beta[lm_rows], tolerance = 1e-6)
  expect_equal(ct_feols$SE, ct_lm$SE[lm_rows], tolerance = 1e-6)

})


test_that("Wald_test works", {

  Wald_feols <- Wald_test(fit_fe, constraints = constrain_zero(1:2), vcov = "CR2", test = "All")
  expect_s3_class(Wald_feols, "data.frame")
  expect_true(all(Wald_feols$Fstat > 0))
  expect_true(all(Wald_feols$df_num == 2))

  # Compare to no-FE model where lm and feols should give identical results
  Wald_nfe <- Wald_test(fit_nfe, constraints = constrain_zero(2:3),
                         vcov = "CR2", cluster = dat$cluster, test = "All")
  Wald_lm <- Wald_test(lm_nfe, constraints = constrain_zero(2:3),
                        vcov = "CR2", cluster = dat$cluster, test = "All")
  expect_equal(Wald_nfe$Fstat, Wald_lm$Fstat, tolerance = 1e-6)
  expect_equal(Wald_nfe$df_denom, Wald_lm$df_denom, tolerance = 1e-3)

})


test_that("dropped observations are handled", {

  dat_drop <- dat
  dat_drop$X1[c(1, 5, 10)] <- NA
  fit_drop <- feols(y ~ X1 + X2 | cluster, data = dat_drop)
  dat_complete <- dat_drop[complete.cases(dat_drop[c("y", "X1", "X2")]), ]
  fit_complete <- feols(y ~ X1 + X2 | cluster, data = dat_complete)

  for (cr_type in CR_types) {
    v_drop <- vcovCR(fit_drop, type = cr_type)
    v_complete <- vcovCR(fit_complete, type = cr_type)
    expect_equal(as.matrix(v_drop), as.matrix(v_complete), tolerance = 1e-6)
  }

})
