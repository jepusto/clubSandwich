context("2-level lme objects")
suppressMessages(library(lme4, quietly=TRUE))
library(nlme, quietly=TRUE, warn.conflicts=FALSE)
library(mlmRev, quietly=TRUE, warn.conflicts=FALSE)

obj_A <- lme(weight ~ Time * Diet, data=BodyWeight, ~ Time | Rat)
obj_A2 <- update(obj_A, weights = varPower())
obj_A3 <- update(obj_A, correlation = corExp(form = ~ Time))
obj_A4 <- update(obj_A2, correlation = corExp(form = ~ Time))
obj_B <- lme(distance ~ age, random = ~ age, data = Orthodont)

obj <- obj_A2

test_that("bread works", {
  expect_true(check_bread(obj_A, cluster = BodyWeight$Rat, y = BodyWeight$weight))
  expect_true(check_bread(obj_A2, cluster = BodyWeight$Rat, y = BodyWeight$weight, tol = 5 * 10^-5))
  expect_true(check_bread(obj_A3, cluster = BodyWeight$Rat, y = BodyWeight$weight))
  expect_true(check_bread(obj_A4, cluster = BodyWeight$Rat, y = BodyWeight$weight))
  expect_true(check_bread(obj_B, cluster = Orthodont$Subject, y = Orthodont$distance))
  
  expect_equal(vcov(obj_A), obj_A$sigma^2 * bread(obj_A) / v_scale(obj_A))
  expect_equal(vcov(obj_A2), obj_A2$sigma^2 * bread(obj_A2) / v_scale(obj_A2))
  expect_equal(vcov(obj_A3), obj_A3$sigma^2 * bread(obj_A3) / v_scale(obj_A3))
  expect_equal(vcov(obj_A4), obj_A4$sigma^2 * bread(obj_A4) / v_scale(obj_A4))
  expect_equal(vcov(obj_B), obj_B$sigma^2 * bread(obj_B) / v_scale(obj_B))

  })

test_that("vcovCR options work for CR2", {
  CR2_A <- vcovCR(obj_A, type = "CR2")
  expect_identical(vcovCR(obj_A, cluster = BodyWeight$Rat, type = "CR2"), CR2_A)
  expect_identical(vcovCR(obj_A, type = "CR2", inverse_var = TRUE), CR2_A)
  expect_false(identical(vcovCR(obj_A, type = "CR2", inverse_var = FALSE), CR2_A))
  
  target <- targetVariance(obj_A)
  expect_equal(vcovCR(obj_A, type = "CR2", target = target, inverse_var = TRUE), CR2_A)
  attr(CR2_A, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A, type = "CR2", target = target, inverse_var = FALSE), CR2_A)
  
  CR2_A2 <- vcovCR(obj_A2, type = "CR2")
  expect_identical(vcovCR(obj_A2, cluster = BodyWeight$Rat, type = "CR2"), CR2_A2)
  expect_identical(vcovCR(obj_A2, type = "CR2", inverse_var = TRUE), CR2_A2)
  expect_false(identical(vcovCR(obj_A2, type = "CR2", inverse_var = FALSE), CR2_A2))
  
  target <- targetVariance(obj_A2)
  expect_equal(vcovCR(obj_A2, type = "CR2", target = target, inverse_var = TRUE), CR2_A2)
  attr(CR2_A2, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A2, type = "CR2", target = target, inverse_var = FALSE), CR2_A2)
  
  CR2_A3 <- vcovCR(obj_A3, type = "CR2")
  expect_identical(vcovCR(obj_A3, cluster = BodyWeight$Rat, type = "CR2"), CR2_A3)
  expect_identical(vcovCR(obj_A3, type = "CR2", inverse_var = TRUE), CR2_A3)
  expect_false(identical(vcovCR(obj_A3, type = "CR2", inverse_var = FALSE), CR2_A3))
  
  target <- targetVariance(obj_A3)
  expect_equal(vcovCR(obj_A3, type = "CR2", target = target, inverse_var = TRUE), CR2_A3)
  attr(CR2_A3, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A3, type = "CR2", target = target, inverse_var = FALSE), CR2_A3)

  CR2_B <- vcovCR(obj_B, type = "CR2")
  expect_identical(vcovCR(obj_B, cluster = Orthodont$Subject, type = "CR2"), CR2_B)
  expect_identical(vcovCR(obj_B, type = "CR2", inverse_var = TRUE), CR2_B)
  expect_false(identical(vcovCR(obj_B, type = "CR2", inverse_var = FALSE), CR2_B))
  
  target <- targetVariance(obj_B)
  expect_equal(vcovCR(obj_B, type = "CR2", target = target, inverse_var = TRUE), CR2_B)
  attr(CR2_B, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_B, type = "CR2", target = target, inverse_var = FALSE), CR2_B)
})

test_that("vcovCR options work for CR4", {
  CR4_A <- vcovCR(obj_A, type = "CR4")
  expect_identical(vcovCR(obj_A, cluster = BodyWeight$Rat, type = "CR4"), CR4_A)
  expect_identical(vcovCR(obj_A, type = "CR4", inverse_var = TRUE), CR4_A)
  expect_false(identical(vcovCR(obj_A, type = "CR4", inverse_var = FALSE), CR4_A))
  
  target <- targetVariance(obj_A)
  expect_equal(vcovCR(obj_A, type = "CR4", target = target, inverse_var = TRUE), CR4_A)
  attr(CR4_A, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A, type = "CR4", target = target, inverse_var = FALSE), CR4_A)
  
  CR4_B <- vcovCR(obj_B, type = "CR4")
  expect_identical(vcovCR(obj_B, cluster = Orthodont$Subject, type = "CR4"), CR4_B)
  expect_identical(vcovCR(obj_B, type = "CR4", inverse_var = TRUE), CR4_B)
  expect_false(identical(vcovCR(obj_B, type = "CR4", inverse_var = FALSE), CR4_B))
  
  target <- targetVariance(obj_B)
  expect_equal(vcovCR(obj_B, type = "CR4", target = target, inverse_var = TRUE), CR4_B)
  attr(CR4_B, "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_B, type = "CR4", target = target, inverse_var = FALSE), CR4_B)
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(obj_A, vcov = "CR2"))
  expect_true(check_CR(obj_B, vcov = "CR2"))
  expect_true(check_CR(obj_A, vcov = "CR4"))
  expect_true(check_CR(obj_B, vcov = "CR4"))
})


CR_types <- paste0("CR",0:4)

test_that("Order doesn't matter.", {
  re_order <- sample(nrow(BodyWeight))
  dat_scramble <- BodyWeight[re_order,]
  obj_scramble <- update(obj_A, data = dat_scramble)
  
  CR_fit <- lapply(CR_types, function(x) vcovCR(obj_A, type = x))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(obj_scramble, type = x))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(obj_A, vcov = x, test = "All"))
  test_scramble <- lapply(CR_types, function(x) coef_test(obj_scramble, vcov = x, test = "All"))
  expect_equal(test_fit, test_scramble, tolerance = 10^-6)
  
  constraints <- combn(length(coef(obj_A)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(obj_A, constraints = constraints, vcov = "CR2", test = "All")
  Wald_scramble <- Wald_test(obj_scramble, constraints = constraints, vcov = "CR2", test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})


test_that("clubSandwich works with dropped observations", {
  dat_miss <- BodyWeight
  dat_miss$weight[sample.int(nrow(BodyWeight), size = round(nrow(BodyWeight) / 10))] <- NA
  obj_dropped <- update(obj_A, data = dat_miss, na.action = na.omit)
  obj_complete <- update(obj_A, data = dat_miss, subset = !is.na(weight))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(obj_dropped, vcov = x, test = "All"))
  test_complete <- lapply(CR_types, function(x) coef_test(obj_complete, vcov = x, test = "All"))
  expect_identical(test_drop, test_complete)
})



test_that("lme agrees with gls", {
  lme_fit <- lme(weight ~ Time * Diet, data=BodyWeight, ~ 1 | Rat)
  gls_fit <- gls(weight ~ Time * Diet, data=BodyWeight, 
                 correlation = corCompSymm(form = ~ 1 | Rat))
  
  CR_lme <- lapply(CR_types, function(x) vcovCR(lme_fit, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_fit, type = x))
  max_ratio <- mapply(function(a, b) max(abs(a / b - 1)), CR_lme, CR_gls)
  expect_true(all(max_ratio < 10^-8))
  
  test_lme <- lapply(CR_types, function(x) coef_test(lme_fit, vcov = x, test = "All"))
  test_gls <- lapply(CR_types, function(x) coef_test(gls_fit, vcov = x, test = "All"))
  compare_tests <- mapply(function(a, b) max(abs(a / b - 1), na.rm = TRUE), test_lme, test_gls)
  expect_true(all(compare_tests < 10^-8))
  
  constraints <- c(combn(length(coef(lme_fit)), 2, simplify = FALSE),
                   combn(length(coef(lme_fit)), 3, simplify = FALSE))
  Wald_lme <- Wald_test(lme_fit, constraints = constraints, vcov = "CR2", test = "All")
  Wald_gls <- Wald_test(gls_fit, constraints = constraints, vcov = "CR2", test = "All")
  expect_equal(Wald_lme, Wald_gls)
})



test_that("Emply levels are dropped in model_matrix", {
  
  data(AchievementAwardsRCT)
  
  AA_RCT_females <- subset(AchievementAwardsRCT, sex=="Girl" & year != "1999")
  
  AA_RCT_females <- within(AA_RCT_females, {
    sibs_4 <- siblings >= 4
    treated2001 <- treated * (year=="2001")
  })
  
  lme_fit <- lme(Bagrut_status ~ year * school_type + 
                   father_ed + mother_ed + immigrant + sibs_4 + 
                   qrtl + treated2001:half, 
                 random = ~ 1 | school_id, 
                 data = AA_RCT_females)
  
  betas <- fixef(lme_fit)
  X <- model_matrix(lme_fit)
  expect_identical(names(betas), colnames(X))
  
})
