context("mmrm objects")
set.seed(20190513)

skip_if_not_installed("mmrm")
skip_if_not_installed("nlme")

library(mmrm, quietly = TRUE, warn.conflicts = FALSE)

# Dataset 1: fev_data (from mmrm package)

data(fev_data, package = "mmrm")

# Unstructured covariance
obj_us <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)

# Toeplitz covariance
obj_toep <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + toep(AVISIT | USUBJID),
  data = fev_data
)

# Grouped covariance (separate covariance per treatment arm)
obj_grouped <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | ARMCD / USUBJID),
  data = fev_data
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

# Compound symmetry covariance
obj_bw <- mmrm(
  weight ~ Diet + cs(Time_f | Rat),
  data = bw_data
)

ff_bw <- component(obj_bw, "full_frame")
cluster_bw <- droplevels(as.factor(ff_bw[[obj_bw$formula_parts$subject_var]]))


# Dataset 3: Orthodont (from nlme) — orthodontic growth data
# 27 subjects (16 Male, 11 Female), ages 8/10/12/14

data(Orthodont, package = "nlme")
orth_data <- as.data.frame(Orthodont)
orth_data$age_f <- factor(orth_data$age)
orth_data$Subject <- factor(as.character(orth_data$Subject)) # Orthodont$Subject is an ordered factor

# AR(1) covariance (interaction breaks perfect orthogonality in design)
obj_orth <- mmrm(
  distance ~ Sex * age_f + ar1(age_f | Subject),
  data = orth_data
)
# Use the Sex * age f interaction (rather than just Sex + age f) to break perfect
# orthogonality in the design matrix, which avoids numerical issues in the bread check

ff_orth <- component(obj_orth, "full_frame")
cluster_orth <- droplevels(as.factor(ff_orth[[obj_orth$formula_parts$subject_var]]))


### Tests

test_that("bread works", {
  # fev_data models
  expect_true(check_bread(obj_us, cluster = cluster_us, y = ff_us$FEV1))
  expect_true(check_bread(obj_toep, cluster = cluster_toep, y = ff_toep$FEV1))

  # BodyWeight model
  expect_true(check_bread(obj_bw, cluster = cluster_bw, y = ff_bw$weight))

  # Orthodont model
  expect_true(check_bread(obj_orth, cluster = cluster_orth, y = ff_orth$distance))

  # Verify bread formula: vcov(obj) = bread(obj) / v_scale(obj)
  expect_equal(vcov(obj_us), bread(obj_us) / v_scale(obj_us))
  expect_equal(vcov(obj_toep), bread(obj_toep) / v_scale(obj_toep))
  expect_equal(vcov(obj_bw), bread(obj_bw) / v_scale(obj_bw))
  expect_equal(vcov(obj_orth), bread(obj_orth) / v_scale(obj_orth))
})


test_that("vcovCR options work for CR2", {
  # fev_data: unstructured
  CR2_us <- vcovCR(obj_us, type = "CR2")
  expect_equal(vcovCR(obj_us, cluster = cluster_us, type = "CR2"), CR2_us)
  expect_equal(vcovCR(obj_us, type = "CR2", inverse_var = TRUE), CR2_us)
  expect_false(identical(vcovCR(obj_us, type = "CR2", inverse_var = FALSE), CR2_us))

  target <- targetVariance(obj_us, cluster_us)
  expect_equal(vcovCR(obj_us, type = "CR2", target = target, inverse_var = TRUE), CR2_us)
  attr(CR2_us, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_us, type = "CR2", target = target, inverse_var = FALSE), CR2_us)

  # BodyWeight: compound symmetry
  CR2_bw <- vcovCR(obj_bw, type = "CR2")
  expect_equal(vcovCR(obj_bw, cluster = cluster_bw, type = "CR2"), CR2_bw)
  expect_equal(vcovCR(obj_bw, type = "CR2", inverse_var = TRUE), CR2_bw)
  expect_false(identical(vcovCR(obj_bw, type = "CR2", inverse_var = FALSE), CR2_bw))

  target_bw <- targetVariance(obj_bw, cluster_bw)
  expect_equal(vcovCR(obj_bw, type = "CR2", target = target_bw, inverse_var = TRUE), CR2_bw)
  attr(CR2_bw, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_bw, type = "CR2", target = target_bw, inverse_var = FALSE), CR2_bw)

  # Orthodont: AR(1)
  CR2_orth <- vcovCR(obj_orth, type = "CR2")
  expect_equal(vcovCR(obj_orth, cluster = cluster_orth, type = "CR2"), CR2_orth)
  expect_equal(vcovCR(obj_orth, type = "CR2", inverse_var = TRUE), CR2_orth)
  expect_false(identical(vcovCR(obj_orth, type = "CR2", inverse_var = FALSE), CR2_orth))

  target_orth <- targetVariance(obj_orth, cluster_orth)
  expect_equal(vcovCR(obj_orth, type = "CR2", target = target_orth, inverse_var = TRUE), CR2_orth)
  attr(CR2_orth, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_orth, type = "CR2", target = target_orth, inverse_var = FALSE), CR2_orth)
})


test_that("vcovCR options work for CR4", {
  # fev_data
  CR4_us <- vcovCR(obj_us, type = "CR4")
  expect_equal(vcovCR(obj_us, cluster = cluster_us, type = "CR4"), CR4_us)
  expect_equal(vcovCR(obj_us, type = "CR4", inverse_var = TRUE), CR4_us)
  expect_false(identical(vcovCR(obj_us, type = "CR4", inverse_var = FALSE), CR4_us))

  target <- targetVariance(obj_us, cluster_us)
  expect_equal(vcovCR(obj_us, type = "CR4", target = target, inverse_var = TRUE), CR4_us)
  attr(CR4_us, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_us, type = "CR4", target = target, inverse_var = FALSE), CR4_us)

  # Orthodont
  CR4_orth <- vcovCR(obj_orth, type = "CR4")
  expect_equal(vcovCR(obj_orth, cluster = cluster_orth, type = "CR4"), CR4_orth)
  expect_equal(vcovCR(obj_orth, type = "CR4", inverse_var = TRUE), CR4_orth)
  expect_false(identical(vcovCR(obj_orth, type = "CR4", inverse_var = FALSE), CR4_orth))

  target_orth <- targetVariance(obj_orth, cluster_orth)
  expect_equal(vcovCR(obj_orth, type = "CR4", target = target_orth, inverse_var = TRUE), CR4_orth)
  attr(CR4_orth, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_orth, type = "CR4", target = target_orth, inverse_var = FALSE), CR4_orth)
})


test_that("CR2 and CR4 are target-unbiased", {
  # fev_data models
  expect_true(check_CR(obj_us, vcov = "CR2"))
  expect_true(check_CR(obj_us, vcov = "CR4"))
  expect_true(check_CR(obj_toep, vcov = "CR2"))
  expect_true(check_CR(obj_toep, vcov = "CR4"))

  # BodyWeight
  expect_true(check_CR(obj_bw, vcov = "CR2"))
  expect_true(check_CR(obj_bw, vcov = "CR4"))

  # Orthodont
  expect_true(check_CR(obj_orth, vcov = "CR2"))
  expect_true(check_CR(obj_orth, vcov = "CR4"))
})


test_that("grouped covariance model works", {
  # bread works for grouped model
  expect_true(check_bread(obj_grouped, cluster = cluster_grouped, y = ff_grouped$FEV1))
  expect_equal(vcov(obj_grouped), bread(obj_grouped) / v_scale(obj_grouped))

  # CR2 produces a valid vcovCR object
  CR2_grouped <- vcovCR(obj_grouped, type = "CR2")
  expect_s3_class(CR2_grouped, "vcovCR")
  expect_equal(vcovCR(obj_grouped, cluster = cluster_grouped, type = "CR2"), CR2_grouped)

  # CR2 and CR4 are target-unbiased for grouped model
  expect_true(check_CR(obj_grouped, vcov = "CR2"))
  expect_true(check_CR(obj_grouped, vcov = "CR4"))

  # coef_test works on grouped model
  test_grouped <- coef_test(obj_grouped, vcov = "CR2")
  expect_s3_class(test_grouped, "data.frame")
  expect_true(all(test_grouped$SE > 0))
})


test_that("Order doesn't matter", {
  # fev_data
  check_sort_order(obj_us, fev_data)

  # BodyWeight
  check_sort_order(obj_bw, bw_data)

  # Orthodont
  check_sort_order(obj_orth, orth_data)
})


test_that("different covariance structures all work", {
  CR_types <- paste0("CR", 0:4)

  # us (fev_data)
  CR_us <- lapply(CR_types, function(x) vcovCR(obj_us, type = x))
  expect_true(all(sapply(CR_us, inherits, "vcovCR")))

  # toep (fev_data)
  CR_toep <- lapply(CR_types, function(x) vcovCR(obj_toep, type = x))
  expect_true(all(sapply(CR_toep, inherits, "vcovCR")))

  # cs (BodyWeight)
  CR_bw <- lapply(CR_types, function(x) vcovCR(obj_bw, type = x))
  expect_true(all(sapply(CR_bw, inherits, "vcovCR")))

  # ar1 (Orthodont)
  CR_orth <- lapply(CR_types, function(x) vcovCR(obj_orth, type = x))
  expect_true(all(sapply(CR_orth, inherits, "vcovCR")))
})


test_that("coef_test and conf_int work", {
  # fev_data
  test_us <- coef_test(obj_us, vcov = "CR2", test = "All")
  expect_s3_class(test_us, "data.frame")
  expect_true(all(test_us$SE > 0))
  expect_true(all(test_us$df_Satt > 0))

  ci_us <- conf_int(obj_us, vcov = "CR2")
  expect_s3_class(ci_us, "data.frame")
  expect_true(all(ci_us$CI_L < ci_us$CI_U))
  expect_true(all(ci_us$SE > 0))

  # BodyWeight
  test_bw <- coef_test(obj_bw, vcov = "CR2", test = "All")
  expect_s3_class(test_bw, "data.frame")
  expect_true(all(test_bw$SE > 0))

  ci_bw <- conf_int(obj_bw, vcov = "CR2")
  expect_s3_class(ci_bw, "data.frame")
  expect_true(all(ci_bw$CI_L < ci_bw$CI_U))

  # Orthodont
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


test_that("Wald_test works", {
  # fev_data — test first 3 pairwise contrasts
  coefs_us <- coef_CS(obj_us)
  pairs_us <- utils::combn(length(coefs_us), 2, simplify = FALSE)[1:3]
  constraints_us <- lapply(pairs_us, constrain_zero, coefs = coefs_us)

  Wald_us <- Wald_test(obj_us, constraints = constraints_us, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_us))
  expect_true(all(sapply(Wald_us, function(w) all(w$Fstat >= 0))))

  # BodyWeight — all pairwise (only 3 coefficients, so 3 pairs)
  coefs_bw <- coef_CS(obj_bw)
  pairs_bw <- utils::combn(length(coefs_bw), 2, simplify = FALSE)
  constraints_bw <- lapply(pairs_bw, constrain_zero, coefs = coefs_bw)

  Wald_bw <- Wald_test(obj_bw, constraints = constraints_bw, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_bw))
  expect_true(all(sapply(Wald_bw, function(w) all(w$Fstat >= 0))))

  # Orthodont — first 3 pairwise contrasts
  coefs_orth <- coef_CS(obj_orth)
  pairs_orth <- utils::combn(length(coefs_orth), 2, simplify = FALSE)[1:3]
  constraints_orth <- lapply(pairs_orth, constrain_zero, coefs = coefs_orth)

  Wald_orth <- Wald_test(obj_orth, constraints = constraints_orth, vcov = "CR2", test = "All")
  expect_true(is.list(Wald_orth))
  expect_true(all(sapply(Wald_orth, function(w) all(w$Fstat >= 0))))
})

# TODO: test dropped observations / missing visits
