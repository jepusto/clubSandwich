context("lmerMod objects")
set.seed(20191217)

suppressMessages(library(lme4, quietly=TRUE))
library(nlme, quietly=TRUE, warn.conflicts=FALSE)

obj_A1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
obj_A2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)

data(Orthodont, package="nlme")

obj_B1 <- lmer(distance ~ age + (1 | Subject), data=Orthodont)
obj_B2 <- lmer(distance ~ age + (age || Subject), data=Orthodont)

data(egsingle, package = "mlmRev")
egsingle <- within(egsingle, {
  size <- (size - mean(size)) / sd(size)
})
  
obj_C1 <- lmer(math ~ year * size + female + black + hispanic + (1 | schoolid) + (1 | childid),
               data = egsingle)
obj_C2 <- lmer(math ~ year * size + female + black + hispanic + (year | schoolid) + (1 | childid),
               data = egsingle, control = lmerControl(check.conv.grad = .makeCC("ignore", tol = 2e-3, relTol = NULL)))

test_that("bread works", {
  
  expect_true(check_bread(obj_A1, cluster = sleepstudy$Subject, y = sleepstudy$Reaction))
  expect_true(check_bread(obj_A2, cluster = sleepstudy$Subject, y = sleepstudy$Reaction))
  expect_true(check_bread(obj_B1, cluster = Orthodont$Subject, y = Orthodont$distance))
  expect_true(check_bread(obj_B2, cluster = Orthodont$Subject, y = Orthodont$distance))
  expect_true(check_bread(obj_C1, cluster = egsingle$schoolid, y = egsingle$math))
  expect_true(check_bread(obj_C2, cluster = egsingle$schoolid, y = egsingle$math))
  
  expect_equal(as.matrix(vcov(obj_A1)), bread(obj_A1) * getME(obj_A1, "sigma")^2 / v_scale(obj_A1))
  expect_equal(as.matrix(vcov(obj_A2)), bread(obj_A2) * getME(obj_A2, "sigma")^2 / v_scale(obj_A2))
  expect_equal(as.matrix(vcov(obj_B1)), bread(obj_B1) * getME(obj_B1, "sigma")^2 / v_scale(obj_B1))
  expect_equal(as.matrix(vcov(obj_B2)), bread(obj_B2) * getME(obj_B2, "sigma")^2 / v_scale(obj_B2))
  expect_equal(as.matrix(vcov(obj_C1)), bread(obj_C1) * getME(obj_C1, "sigma")^2 / v_scale(obj_C1))
  expect_equal(as.matrix(vcov(obj_C2)), bread(obj_C2) * getME(obj_C2, "sigma")^2 / v_scale(obj_C2))

})

test_that("vcovCR options work for CR2", {
  
  CR2_A <- vcovCR(obj_A1, type = "CR2")
  expect_equal(vcovCR(obj_A1, cluster = sleepstudy$Subject, type = "CR2"), CR2_A)
  expect_equal(vcovCR(obj_A1, type = "CR2", inverse_var = TRUE), CR2_A)
  expect_false(identical(vcovCR(obj_A1, type = "CR2", inverse_var = FALSE), CR2_A))
  target <- targetVariance(obj_A1)
  expect_equal(vcovCR(obj_A1, type = "CR2", target = target, inverse_var = TRUE), CR2_A, check.attributes = FALSE)
  expect_equal(vcovCR(obj_A1, type = "CR2", target = target, inverse_var = FALSE), CR2_A, check.attributes = FALSE)
  
  CR2_B <- vcovCR(obj_B1, type = "CR2")
  expect_equal(vcovCR(obj_B1, cluster = Orthodont$Subject, type = "CR2"), CR2_B)
  expect_equal(vcovCR(obj_B1, type = "CR2", inverse_var = TRUE), CR2_B)
  expect_false(identical(vcovCR(obj_B1, type = "CR2", inverse_var = FALSE), CR2_B))
  target <- targetVariance(obj_B1)
  expect_equal(vcovCR(obj_B1, type = "CR2", target = target, inverse_var = TRUE), CR2_B, check.attributes = FALSE)
  expect_equal(vcovCR(obj_B1, type = "CR2", target = target, inverse_var = FALSE), CR2_B, check.attributes = FALSE)
  
  CR2_C <- vcovCR(obj_C1, type = "CR2")
  expect_equal(vcovCR(obj_C1, cluster = egsingle$schoolid, type = "CR2"), CR2_C)
  expect_equal(vcovCR(obj_C1, type = "CR2", inverse_var = TRUE), CR2_C)
  expect_false(identical(vcovCR(obj_C1, type = "CR2", inverse_var = FALSE), CR2_C))
  target <- targetVariance(obj_C1)
  expect_equal(vcovCR(obj_C1, type = "CR2", target = target, inverse_var = TRUE), CR2_C, check.attributes = FALSE)
  expect_equal(vcovCR(obj_C1, type = "CR2", target = target, inverse_var = FALSE), CR2_C, check.attributes = FALSE)
  
})

test_that("vcovCR options work for CR4", {
  
  CR4_A <- vcovCR(obj_A1, type = "CR4")
  expect_identical(vcovCR(obj_A1, cluster = sleepstudy$Subject, type = "CR4"), CR4_A)
  expect_identical(vcovCR(obj_A1, type = "CR4", inverse_var = TRUE), CR4_A)
  expect_false(identical(vcovCR(obj_A1, type = "CR4", inverse_var = FALSE), CR4_A))
  target <- targetVariance(obj_A1)
  expect_equal(vcovCR(obj_A1, type = "CR4", target = target, inverse_var = TRUE), CR4_A, check.attributes = FALSE)
  expect_equal(vcovCR(obj_A1, type = "CR4", target = target, inverse_var = FALSE), CR4_A, check.attributes = FALSE)
  
  CR4_B <- vcovCR(obj_B1, type = "CR4")
  expect_identical(vcovCR(obj_B1, cluster = Orthodont$Subject, type = "CR4"), CR4_B)
  expect_identical(vcovCR(obj_B1, type = "CR4", inverse_var = TRUE), CR4_B)
  expect_false(identical(vcovCR(obj_B1, type = "CR4", inverse_var = FALSE), CR4_B))
  target <- targetVariance(obj_B1)
  expect_equal(vcovCR(obj_B1, type = "CR4", target = target, inverse_var = TRUE), CR4_B, check.attributes = FALSE)
  expect_equal(vcovCR(obj_B1, type = "CR4", target = target, inverse_var = FALSE), CR4_B, check.attributes = FALSE)
  
  CR4_C <- vcovCR(obj_C1, type = "CR4")
  expect_identical(vcovCR(obj_C1, cluster = egsingle$schoolid, type = "CR4"), CR4_C)
  expect_identical(vcovCR(obj_C1, type = "CR4", inverse_var = TRUE), CR4_C)
  expect_false(identical(vcovCR(obj_C1, type = "CR4", inverse_var = FALSE), CR4_C))
  target <- targetVariance(obj_C1)
  expect_equal(vcovCR(obj_C1, type = "CR4", target = target, inverse_var = TRUE), CR4_C, check.attributes = FALSE)
  expect_equal(vcovCR(obj_C1, type = "CR4", target = target, inverse_var = FALSE), CR4_C, check.attributes = FALSE)
  
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(obj_A1, vcov = "CR2"))
  expect_true(check_CR(obj_A2, vcov = "CR2"))
  expect_true(check_CR(obj_B1, vcov = "CR2"))
  expect_true(check_CR(obj_B2, vcov = "CR2"))
  expect_true(check_CR(obj_C1, vcov = "CR2"))
  expect_true(check_CR(obj_C2, vcov = "CR2"))
  expect_true(check_CR(obj_A1, vcov = "CR4"))
  expect_true(check_CR(obj_A2, vcov = "CR4"))
  expect_true(check_CR(obj_B1, vcov = "CR4"))
  expect_true(check_CR(obj_B2, vcov = "CR4"))
  expect_true(check_CR(obj_C1, vcov = "CR4"))
  expect_true(check_CR(obj_C2, vcov = "CR4"))
})


CR_types <- paste0("CR",0:4)

test_that("Order doesn't matter.", {
  
  skip_on_cran()
  
  # Model A1
  check_sort_order(obj = obj_A1, dat = sleepstudy)

  # Model C1
  check_sort_order(obj = obj_C1, dat = egsingle)
})


test_that("clubSandwich works with dropped observations", {
  
  dat_miss <- sleepstudy
  dat_miss$Reaction[sample.int(nrow(sleepstudy), size = round(nrow(sleepstudy) / 10))] <- NA
  obj_dropped <- update(obj_A1, data = dat_miss, na.action = na.omit)
  obj_complete <- update(obj_A1, data = dat_miss, subset = !is.na(Reaction))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(obj_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(obj_complete, vcov = x, test = "All", p_values = FALSE))
  expect_identical(test_drop, test_complete)
})



test_that("lmer agrees with lme", {
  
  data(BodyWeight, package="nlme")
  
  lmer_fit <- lmer(weight ~ Time * Diet + (1 | Rat), data=BodyWeight)
  lme_fit <- lme(weight ~ Time * Diet, data=BodyWeight, ~ 1 | Rat)
  
  expect_equal(coef_CS(lmer_fit), coef_CS(lme_fit))
  expect_equal(nobs(lmer_fit), nobs(lme_fit))
  expect_equal(model_matrix(lmer_fit), model_matrix(lme_fit), check.attributes = FALSE)
  expect_equal(residuals_CS(lmer_fit), residuals_CS(lme_fit), check.attributes = FALSE)
  expect_equal(v_scale(lmer_fit), v_scale(lme_fit))
  
  p <- length(coef_CS(lmer_fit))
  expect_equal(bread(lmer_fit) / bread(lme_fit), matrix(1, p, p), check.attributes = FALSE, tol = 10^-6)
  expect_equal(targetVariance(lmer_fit), targetVariance(lme_fit), check.attributes = FALSE, tol = 10^-6)
  expect_equal(weightMatrix(lmer_fit), weightMatrix(lme_fit), check.attributes = FALSE, tol = 10^-6)
  
  CR_lmer <- lapply(CR_types, function(x) vcovCR(lmer_fit, type = x))
  CR_lme <- lapply(CR_types, function(x) vcovCR(lme_fit, type = x))
  expect_equivalent(CR_lmer, CR_lme, tolerance = 10^-6)
  
  test_lmer <- lapply(CR_types, function(x) coef_test(lmer_fit, vcov = x, test = "All", p_values = FALSE))
  test_lme <- lapply(CR_types, function(x) coef_test(lme_fit, vcov = x, test = "All", p_values = FALSE))
  compare_ttests(test_lmer, test_lme, tol = 10^-10)
  
  constraints <- c(combn(length(coef_CS(lmer_fit)), 2, simplify = FALSE),
                   combn(length(coef_CS(lmer_fit)), 3, simplify = FALSE))
  Wald_lmer <- Wald_test(lmer_fit, constraints = constrain_zero(constraints), vcov = "CR2", test = "All")
  Wald_lme <- Wald_test(lme_fit, constraints = constrain_zero(constraints), vcov = "CR2", test = "All")
  compare_Waldtests(Wald_lmer, Wald_lme)
  
})


test_that("Emply levels are dropped in model_matrix", {
  
  data(AchievementAwardsRCT)
  
  AA_RCT_females <- subset(AchievementAwardsRCT, sex=="Girl" & year != "1999")
  
  AA_RCT_females <- within(AA_RCT_females, {
    sibs_4 <- siblings >= 4
    treated2001 <- treated * (year=="2001")
  })
  
  lmer_fit <- lmer(Bagrut_status ~ year * school_type + 
                     father_ed + mother_ed + immigrant + sibs_4 + 
                     qrtl + treated2001:half + (1 | school_id),
                   data = AA_RCT_females)
  
  betas <- fixef(lmer_fit)
  X <- model_matrix(lmer_fit)
  expect_identical(names(betas), colnames(X))
  
})



test_that("Possible to cluster at higher level than random effects", {
  
  n_districts <- 10
  n_schools_per <- rnbinom(n_districts, size = 4, prob = 0.3)
  n_schools <- sum(n_schools_per)
  n_students_per <- 10
  n_students <- n_schools * n_students_per
  
  # identifiers for each level
  district_id <- factor(rep(1:n_districts, n_schools_per * n_students_per))
  school_id <- factor(rep(1:sum(n_schools_per), each = n_students_per))
  student_id <- 1:n_students
  
  # simulated outcome
  Y <- rnorm(n_districts)[district_id] + rnorm(n_schools)[school_id] + rnorm(n_students)
  X <- rnorm(n_students)
  dat <- data.frame(district_id, school_id, student_id, Y, X)
  dat_scramble <- dat[sample(nrow(dat)),]
  
  # fit two-level model
  lme_2level <- lmer(Y ~ X + (1 | school_id), data = dat)
  
  # cluster at level 3
  V <- vcovCR(lme_2level, type = "CR2", cluster = dat$district_id)
  expect_is(V, "vcovCR")
  expect_error(vcovCR(lme_2level, type = "CR2", cluster = dat_scramble$district_id))
  
  # check that result does not depend on sort-order
  V_scramble <- vcovCR(lmer(Y ~ X + (1 | school_id), data = dat_scramble), 
                       type = "CR2", cluster = dat_scramble$district_id)
  expect_equal(as.matrix(V), as.matrix(V_scramble))
})
