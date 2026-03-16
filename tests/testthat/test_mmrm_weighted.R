context("weighted mmrm objects")
set.seed(20190513)

skip_if_not_installed("mmrm")
skip_if_not_installed("nlme")

library(mmrm, quietly = TRUE, warn.conflicts = FALSE)

# Dataset 1: fev_data (from mmrm package)

data(fev_data, package = "mmrm")

fev_data$wt <- 1L + rpois(nrow(fev_data), lambda = 3)

# Unstructured covariance with weights
obj_us <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
  data = fev_data,
  weights = fev_data$wt
)

# Toeplitz covariance with weights
obj_toep <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + toep(AVISIT | USUBJID),
  data = fev_data,
  weights = fev_data$wt
)

# Grouped covariance (separate covariance per treatment arm) with weights
obj_grouped <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | ARMCD / USUBJID),
  data = fev_data,
  weights = fev_data$wt
)

# Extract cluster variables for fev_data models
ff_us <- component(obj_us, "full_frame")
cluster_us <- droplevels(as.factor(ff_us[[obj_us$formula_parts$subject_var]]))

ff_toep <- component(obj_toep, "full_frame")
cluster_toep <- droplevels(as.factor(ff_toep[[obj_toep$formula_parts$subject_var]]))

ff_grouped <- component(obj_grouped, "full_frame")
cluster_grouped <- droplevels(as.factor(ff_grouped[[obj_grouped$formula_parts$subject_var]]))


# Dataset 2: BodyWeight (from nlme) — rat growth data
# 16 rats, 3 diet groups, 11 time points

data(BodyWeight, package = "nlme")
bw_data <- as.data.frame(BodyWeight)
bw_data$Time_f <- factor(bw_data$Time)
bw_data$wt <- 1L + rpois(nrow(bw_data), lambda = 2)

# Compound symmetry covariance with weights
obj_bw <- mmrm(
  weight ~ Diet + cs(Time_f | Rat),
  data = bw_data,
  weights = bw_data$wt
)

ff_bw <- component(obj_bw, "full_frame")
cluster_bw <- droplevels(as.factor(ff_bw[[obj_bw$formula_parts$subject_var]]))


# Dataset 3: Orthodont (from nlme) — orthodontic growth data
# 27 subjects (16 Male, 11 Female), ages 8/10/12/14

data(Orthodont, package = "nlme")
orth_data <- as.data.frame(Orthodont)
orth_data$age_f <- factor(orth_data$age)
orth_data$Subject <- factor(as.character(orth_data$Subject))
orth_data$wt <- 1L + rpois(nrow(orth_data), lambda = 2)

# AR(1) covariance (interaction breaks perfect orthogonality in design) with weights
obj_orth <- mmrm(
  distance ~ Sex * age_f + ar1(age_f | Subject),
  data = orth_data,
  weights = orth_data$wt
)

ff_orth <- component(obj_orth, "full_frame")
cluster_orth <- droplevels(as.factor(ff_orth[[obj_orth$formula_parts$subject_var]]))


### Tests

test_that("bread works for weighted models", {
  # fev_data weighted models
  expect_true(check_bread(obj_us, cluster = cluster_us, y = ff_us$FEV1))
  expect_true(check_bread(obj_toep, cluster = cluster_toep, y = ff_toep$FEV1))

  # BodyWeight weighted model
  expect_true(check_bread(obj_bw, cluster = cluster_bw, y = ff_bw$weight))

  # Orthodont weighted model
  expect_true(check_bread(obj_orth, cluster = cluster_orth, y = ff_orth$distance))

  # Verify bread formula: vcov(obj) = bread(obj) / v_scale(obj)
  expect_equal(vcov(obj_us), bread(obj_us) / v_scale(obj_us))
  expect_equal(vcov(obj_toep), bread(obj_toep) / v_scale(obj_toep))
  expect_equal(vcov(obj_bw), bread(obj_bw) / v_scale(obj_bw))
  expect_equal(vcov(obj_orth), bread(obj_orth) / v_scale(obj_orth))
})


test_that("targetVariance incorporates weights correctly", {
  # For the weighted model, targetVariance should differ from the unweighted
  # VarCorr by the weight scaling factor G_i^{-1/2} Sigma G_i^{-1/2}

  tv_us <- targetVariance(obj_us, cluster_us)
  vc_us <- VarCorr(obj_us)
  wts_us <- ff_us[["(weights)"]]
  visits_us <- ff_us[[obj_us$formula_parts$visit_var]]
  visit_levels_us <- levels(visits_us)

  obs_by_subject_us <- split(seq_along(cluster_us), cluster_us)

  # Check each subject's target variance matches manual weight computation
  for (subj in names(obs_by_subject_us)) {
    rows <- obs_by_subject_us[[subj]]
    idx <- match(as.character(visits_us[rows]), visit_levels_us)
    Sigma_subj <- unname(vc_us[idx, idx, drop = FALSE])
    w_inv_sqrt <- 1 / sqrt(wts_us[rows])
    expected <- tcrossprod(w_inv_sqrt) * Sigma_subj
    expect_equal(tv_us[[subj]], expected, tolerance = 1e-10)
  }

  # For an unweighted model, targetVariance should equal VarCorr subsets directly
  obj_bw_unwt <- mmrm(
    weight ~ Diet + cs(Time_f | Rat),
    data = bw_data
  )
  ff_bw_unwt <- component(obj_bw_unwt, "full_frame")
  cluster_bw_unwt <- droplevels(as.factor(ff_bw_unwt[[obj_bw_unwt$formula_parts$subject_var]]))

  tv_bw <- targetVariance(obj_bw_unwt, cluster_bw_unwt)
  vc_bw <- VarCorr(obj_bw_unwt)
  visits_bw <- ff_bw_unwt[[obj_bw_unwt$formula_parts$visit_var]]
  visit_levels_bw <- levels(visits_bw)
  obs_by_subject_bw <- split(seq_along(cluster_bw_unwt), cluster_bw_unwt)

  for (subj in names(obs_by_subject_bw)) {
    rows <- obs_by_subject_bw[[subj]]
    idx <- match(as.character(visits_bw[rows]), visit_levels_bw)
    expected <- unname(vc_bw[idx, idx, drop = FALSE])
    expect_equal(tv_bw[[subj]], expected, tolerance = 1e-10)
  }
})


test_that("vcovCR options work for CR2 with weighted models", {
  # fev_data: weighted unstructured
  CR2_us <- vcovCR(obj_us, type = "CR2")
  expect_equal(vcovCR(obj_us, cluster = cluster_us, type = "CR2"), CR2_us)
  expect_equal(vcovCR(obj_us, type = "CR2", inverse_var = TRUE), CR2_us)
  expect_false(identical(vcovCR(obj_us, type = "CR2", inverse_var = FALSE), CR2_us))

  target <- targetVariance(obj_us, cluster_us)
  expect_equal(vcovCR(obj_us, type = "CR2", target = target, inverse_var = TRUE), CR2_us)
  attr(CR2_us, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_us, type = "CR2", target = target, inverse_var = FALSE), CR2_us)

  # BodyWeight weighted
  CR2_bw <- vcovCR(obj_bw, type = "CR2")
  expect_equal(vcovCR(obj_bw, cluster = cluster_bw, type = "CR2"), CR2_bw)
  expect_equal(vcovCR(obj_bw, type = "CR2", inverse_var = TRUE), CR2_bw)
  expect_false(identical(vcovCR(obj_bw, type = "CR2", inverse_var = FALSE), CR2_bw))

  target_bw <- targetVariance(obj_bw, cluster_bw)
  expect_equal(vcovCR(obj_bw, type = "CR2", target = target_bw, inverse_var = TRUE), CR2_bw)
  attr(CR2_bw, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_bw, type = "CR2", target = target_bw, inverse_var = FALSE), CR2_bw)

  # Orthodont weighted
  CR2_orth <- vcovCR(obj_orth, type = "CR2")
  expect_equal(vcovCR(obj_orth, cluster = cluster_orth, type = "CR2"), CR2_orth)
  expect_equal(vcovCR(obj_orth, type = "CR2", inverse_var = TRUE), CR2_orth)
  expect_false(identical(vcovCR(obj_orth, type = "CR2", inverse_var = FALSE), CR2_orth))

  target_orth <- targetVariance(obj_orth, cluster_orth)
  expect_equal(vcovCR(obj_orth, type = "CR2", target = target_orth, inverse_var = TRUE), CR2_orth)
  attr(CR2_orth, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_orth, type = "CR2", target = target_orth, inverse_var = FALSE), CR2_orth)
})


test_that("vcovCR options work for CR4 with weighted models", {
  # fev_data weighted
  CR4_us <- vcovCR(obj_us, type = "CR4")
  expect_equal(vcovCR(obj_us, cluster = cluster_us, type = "CR4"), CR4_us)
  expect_equal(vcovCR(obj_us, type = "CR4", inverse_var = TRUE), CR4_us)
  expect_false(identical(vcovCR(obj_us, type = "CR4", inverse_var = FALSE), CR4_us))

  target <- targetVariance(obj_us, cluster_us)
  expect_equal(vcovCR(obj_us, type = "CR4", target = target, inverse_var = TRUE), CR4_us)
  attr(CR4_us, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_us, type = "CR4", target = target, inverse_var = FALSE), CR4_us)

  # BodyWeight weighted
  CR4_bw <- vcovCR(obj_bw, type = "CR4")
  expect_equal(vcovCR(obj_bw, cluster = cluster_bw, type = "CR4"), CR4_bw)
  expect_equal(vcovCR(obj_bw, type = "CR4", inverse_var = TRUE), CR4_bw)
  expect_false(identical(vcovCR(obj_bw, type = "CR4", inverse_var = FALSE), CR4_bw))

  target_bw <- targetVariance(obj_bw, cluster_bw)
  expect_equal(vcovCR(obj_bw, type = "CR4", target = target_bw, inverse_var = TRUE), CR4_bw)
  attr(CR4_bw, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_bw, type = "CR4", target = target_bw, inverse_var = FALSE), CR4_bw)

  # Orthodont weighted
  CR4_orth <- vcovCR(obj_orth, type = "CR4")
  expect_equal(vcovCR(obj_orth, cluster = cluster_orth, type = "CR4"), CR4_orth)
  expect_equal(vcovCR(obj_orth, type = "CR4", inverse_var = TRUE), CR4_orth)
  expect_false(identical(vcovCR(obj_orth, type = "CR4", inverse_var = FALSE), CR4_orth))

  target_orth <- targetVariance(obj_orth, cluster_orth)
  expect_equal(vcovCR(obj_orth, type = "CR4", target = target_orth, inverse_var = TRUE), CR4_orth)
  attr(CR4_orth, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_orth, type = "CR4", target = target_orth, inverse_var = FALSE), CR4_orth)
})


test_that("CR2 and CR4 are target-unbiased for weighted models", {
  # fev_data weighted
  expect_true(check_CR(obj_us, vcov = "CR2"))
  expect_true(check_CR(obj_us, vcov = "CR4"))
  expect_true(check_CR(obj_toep, vcov = "CR2"))
  expect_true(check_CR(obj_toep, vcov = "CR4"))

  # BodyWeight weighted
  expect_true(check_CR(obj_bw, vcov = "CR2"))
  expect_true(check_CR(obj_bw, vcov = "CR4"))

  # Orthodont weighted
  expect_true(check_CR(obj_orth, vcov = "CR2"))
  expect_true(check_CR(obj_orth, vcov = "CR4"))
})


test_that("grouped covariance model works for weighted models", {
  # bread works for weighted grouped model
  expect_true(check_bread(obj_grouped, cluster = cluster_grouped, y = ff_grouped$FEV1))
  expect_equal(vcov(obj_grouped), bread(obj_grouped) / v_scale(obj_grouped))

  # CR2 produces a valid vcovCR object
  CR2_grouped <- vcovCR(obj_grouped, type = "CR2")
  expect_s3_class(CR2_grouped, "vcovCR")
  expect_equal(vcovCR(obj_grouped, cluster = cluster_grouped, type = "CR2"), CR2_grouped)

  # CR2 and CR4 are target-unbiased for weighted grouped model
  expect_true(check_CR(obj_grouped, vcov = "CR2"))
  expect_true(check_CR(obj_grouped, vcov = "CR4"))

  # coef_test works on weighted grouped model
  test_grouped <- coef_test(obj_grouped, vcov = "CR2")
  expect_s3_class(test_grouped, "data.frame")
  expect_true(all(test_grouped$SE > 0))
})


test_that("Order doesn't matter for weighted models", {
  # Manual sort-order check for weighted mmrm models.
  # check_sort_order uses update(obj, data = ...) which replays the original
  # weights= argument literally, so we refit manually with scrambled data

  # fev_data weighted
  re_order <- sample(nrow(fev_data))
  fev_scramble <- fev_data[re_order, ]
  obj_us_scramble <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = fev_scramble,
    weights = fev_scramble$wt
  )
  CR2_orig <- vcovCR(obj_us, type = "CR2")
  CR2_scramble <- vcovCR(obj_us_scramble, type = "CR2")
  expect_equivalent(as.matrix(CR2_orig), as.matrix(CR2_scramble), tolerance = 1e-6)

  # BodyWeight weighted
  re_order_bw <- sample(nrow(bw_data))
  bw_scramble <- bw_data[re_order_bw, ]
  obj_bw_scramble <- mmrm(
    weight ~ Diet + cs(Time_f | Rat),
    data = bw_scramble,
    weights = bw_scramble$wt
  )
  CR2_bw_orig <- vcovCR(obj_bw, type = "CR2")
  CR2_bw_scramble <- vcovCR(obj_bw_scramble, type = "CR2")
  expect_equivalent(as.matrix(CR2_bw_orig), as.matrix(CR2_bw_scramble), tolerance = 1e-6)

  # Orthodont weighted
  re_order_orth <- sample(nrow(orth_data))
  orth_scramble <- orth_data[re_order_orth, ]
  obj_orth_scramble <- mmrm(
    distance ~ Sex * age_f + ar1(age_f | Subject),
    data = orth_scramble,
    weights = orth_scramble$wt
  )
  CR2_orth_orig <- vcovCR(obj_orth, type = "CR2")
  CR2_orth_scramble <- vcovCR(obj_orth_scramble, type = "CR2")
  expect_equivalent(as.matrix(CR2_orth_orig), as.matrix(CR2_orth_scramble), tolerance = 1e-6)
})


test_that("different covariance structures all work for weighted models", {
  CR_types <- paste0("CR", 0:4)

  # us weighted (fev_data)
  CR_us <- lapply(CR_types, function(x) vcovCR(obj_us, type = x))
  expect_true(all(sapply(CR_us, inherits, "vcovCR")))

  # toep weighted (fev_data)
  CR_toep <- lapply(CR_types, function(x) vcovCR(obj_toep, type = x))
  expect_true(all(sapply(CR_toep, inherits, "vcovCR")))

  # cs weighted (BodyWeight)
  CR_bw <- lapply(CR_types, function(x) vcovCR(obj_bw, type = x))
  expect_true(all(sapply(CR_bw, inherits, "vcovCR")))

  # ar1 weighted (Orthodont)
  CR_orth <- lapply(CR_types, function(x) vcovCR(obj_orth, type = x))
  expect_true(all(sapply(CR_orth, inherits, "vcovCR")))
})


test_that("coef_test and conf_int work for weighted models", {
  # fev_data weighted
  test_us <- coef_test(obj_us, vcov = "CR2", test = "All")
  expect_s3_class(test_us, "data.frame")
  expect_true(all(test_us$SE > 0))
  expect_true(all(test_us$df_Satt > 0))

  ci_us <- conf_int(obj_us, vcov = "CR2")
  expect_s3_class(ci_us, "data.frame")
  expect_true(all(ci_us$CI_L < ci_us$CI_U))
  expect_true(all(ci_us$SE > 0))

  # BodyWeight weighted
  test_bw <- coef_test(obj_bw, vcov = "CR2", test = "All")
  expect_s3_class(test_bw, "data.frame")
  expect_true(all(test_bw$SE > 0))

  ci_bw <- conf_int(obj_bw, vcov = "CR2")
  expect_s3_class(ci_bw, "data.frame")
  expect_true(all(ci_bw$CI_L < ci_bw$CI_U))

  # Orthodont weighted
  test_orth <- coef_test(obj_orth, vcov = "CR2", test = "All")
  expect_s3_class(test_orth, "data.frame")
  expect_true(all(test_orth$SE > 0))

  ci_orth <- conf_int(obj_orth, vcov = "CR2")
  expect_s3_class(ci_orth, "data.frame")
  expect_true(all(ci_orth$CI_L < ci_orth$CI_U))

  # Multiple CR types all produce valid coef_test / conf_int results
  for (cr_type in c("CR0", "CR1", "CR2")) {
    test_cr <- coef_test(obj_us, vcov = cr_type)
    expect_s3_class(test_cr, "data.frame")
    ci_cr <- conf_int(obj_us, vcov = cr_type)
    expect_s3_class(ci_cr, "data.frame")
  }
})


test_that("Wald_test works for weighted models", {
  # fev_data weighted
  coefs_us <- coef_CS(obj_us)
  pairs_us <- utils::combn(length(coefs_us), 2, simplify = FALSE)[1:3]
  constraints_us <- lapply(pairs_us, constrain_zero, coefs = coefs_us)

  Wald_us <- Wald_test(obj_us, constraints = constraints_us, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_us))
  expect_true(all(sapply(Wald_us, function(w) all(w$Fstat >= 0))))

  # BodyWeight weighted
  coefs_bw <- coef_CS(obj_bw)
  pairs_bw <- utils::combn(length(coefs_bw), 2, simplify = FALSE)
  constraints_bw <- lapply(pairs_bw, constrain_zero, coefs = coefs_bw)

  Wald_bw <- Wald_test(obj_bw, constraints = constraints_bw, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_bw))
  expect_true(all(sapply(Wald_bw, function(w) all(w$Fstat >= 0))))

  # Orthodont weighted
  coefs_orth <- coef_CS(obj_orth)
  pairs_orth <- utils::combn(length(coefs_orth), 2, simplify = FALSE)[1:3]
  constraints_orth <- lapply(pairs_orth, constrain_zero, coefs = coefs_orth)

  Wald_orth <- Wald_test(obj_orth, constraints = constraints_orth, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_orth))
  expect_true(all(sapply(Wald_orth, function(w) all(w$Fstat >= 0))))
})


test_that("weighted and unweighted results differ", {
  # Fit identical model formula, one with weights and one without
  obj_us_unwt <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = fev_data
  )

  CR2_wt <- vcovCR(obj_us, type = "CR2")
  CR2_unwt <- vcovCR(obj_us_unwt, type = "CR2")

  # The weighted and unweighted CRVEs should not be identical
  expect_false(identical(as.matrix(CR2_wt), as.matrix(CR2_unwt)))

  # But both should be valid vcovCR objects
  expect_s3_class(CR2_wt, "vcovCR")
  expect_s3_class(CR2_unwt, "vcovCR")
})


test_that("weight scale does not matter for weighted models", {
  # Fit model with scaled-up weights (multiply all weights by a constant)
  fev_data$wt_scaled <- fev_data$wt * 5

  obj_us_scaled <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = fev_data,
    weights = fev_data$wt_scaled
  )

  CR_types <- paste0("CR", 0:4)

  CR_orig <- lapply(CR_types, function(x) vcovCR(obj_us, type = x))
  CR_scaled <- lapply(CR_types, function(x) vcovCR(obj_us_scaled, type = x))

  expect_equal(
    lapply(CR_orig, as.matrix),
    lapply(CR_scaled, as.matrix),
    tolerance = 5e-4
  )
})


test_that("clubSandwich works with dropped observations (weighted)", {

  dat_miss <- fev_data
  dat_miss$FEV1[sample.int(nrow(fev_data), size = round(nrow(fev_data) / 10))] <- NA

  # mmrm auto-drops NA rows (no na.action parameter)
  obj_dropped <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = dat_miss,
    weights = dat_miss$wt
  )

  # mmrm has no subset parameter — manually subset the data
  dat_complete <- dat_miss[!is.na(dat_miss$FEV1), ]
  obj_complete <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = dat_complete,
    weights = dat_complete$wt
  )

  CR_types <- paste0("CR", 0:4)

  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_equal(CR_drop, CR_complete)

  test_drop <- lapply(CR_types, function(x) coef_test(obj_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(obj_complete, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})
