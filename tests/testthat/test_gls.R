context("gls objects")

library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Ovary, package = "nlme")

Ovary$time_int <- 1:nrow(Ovary)

lm_hom <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary)
lm_power <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
                weights = varPower())
lm_AR1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary,
              correlation = corAR1(form = ~ time_int | Mare))
lm_AR1_power <- update(lm_AR1, weights = varPower())

obj <- lm_hom
V <- targetVariance(obj, cluster = Ovary$Mare)

test_that("bread works", {
  expect_true(check_bread(lm_hom, cluster = Ovary$Mare, y = Ovary$follicles))
  expect_true(check_bread(lm_power, cluster = Ovary$Mare, y = Ovary$follicles))
  expect_true(check_bread(lm_AR1, cluster = Ovary$Mare, y = Ovary$follicles))
  expect_true(check_bread(lm_AR1_power, cluster = Ovary$Mare, y = Ovary$follicles))
  
  expect_equal(vcov(lm_hom), lm_hom$sigma^2 * bread(lm_hom) / v_scale(lm_hom))
  expect_equal(vcov(lm_power), lm_power$sigma^2 * bread(lm_power) / v_scale(lm_power))
  expect_equal(vcov(lm_AR1), lm_AR1$sigma^2 * bread(lm_AR1) / v_scale(lm_AR1))
  expect_equal(vcov(lm_AR1_power), lm_AR1_power$sigma^2 * bread(lm_AR1_power) / v_scale(lm_AR1_power))
})

test_that("vcovCR options work for CR2", {
  CR2_AR1 <- vcovCR(lm_AR1, type = "CR2")
  expect_identical(vcovCR(lm_AR1, cluster = Ovary$Mare, type = "CR2"), CR2_AR1)
  expect_identical(vcovCR(lm_AR1, type = "CR2", inverse_var = TRUE), CR2_AR1)
  expect_false(identical(vcovCR(lm_AR1, type = "CR2", inverse_var = FALSE), CR2_AR1))
  
  target <- targetVariance(lm_AR1)
  expect_equal(vcovCR(lm_AR1, type = "CR2", target = target, inverse_var = TRUE), CR2_AR1)
  attr(CR2_AR1, "inverse_var") <- FALSE
  expect_equal(vcovCR(lm_AR1, type = "CR2", target = target, inverse_var = FALSE), CR2_AR1)

  CR2_power <- vcovCR(lm_AR1_power, type = "CR2")
  expect_identical(vcovCR(lm_AR1_power, cluster = Ovary$Mare, type = "CR2"), CR2_power)
  expect_identical(vcovCR(lm_AR1_power, type = "CR2", inverse_var = TRUE), CR2_power)
  expect_false(identical(vcovCR(lm_AR1_power, type = "CR2", inverse_var = FALSE), CR2_power))
  
  target <- targetVariance(lm_AR1_power, cluster = Ovary$Mare)
  expect_equal(vcovCR(lm_AR1_power, type = "CR2", target = target, inverse_var = TRUE), CR2_power)
  attr(CR2_power, "inverse_var") <- FALSE
  expect_equal(vcovCR(lm_AR1_power, type = "CR2", target = target, inverse_var = FALSE), CR2_power)
})


test_that("vcovCR options work for CR4", {
  CR4_AR1 <- vcovCR(lm_AR1, type = "CR4")
  expect_identical(vcovCR(lm_AR1, cluster = Ovary$Mare, type = "CR4"), CR4_AR1)
  expect_identical(vcovCR(lm_AR1, type = "CR4", inverse_var = TRUE), CR4_AR1)
  expect_false(identical(vcovCR(lm_AR1, type = "CR4", inverse_var = FALSE), CR4_AR1))
  
  target <- targetVariance(lm_AR1)
  expect_equal(vcovCR(lm_AR1, type = "CR4", target = target, inverse_var = TRUE), CR4_AR1)
  attr(CR4_AR1, "inverse_var") <- FALSE
  expect_equal(vcovCR(lm_AR1, type = "CR4", target = target, inverse_var = FALSE), CR4_AR1)
  
  CR4_power <- vcovCR(lm_AR1_power, type = "CR4")
  expect_identical(vcovCR(lm_AR1_power, cluster = Ovary$Mare, type = "CR4"), CR4_power)
  expect_identical(vcovCR(lm_AR1_power, type = "CR4", inverse_var = TRUE), CR4_power)
  expect_false(identical(vcovCR(lm_AR1_power, type = "CR4", inverse_var = FALSE), CR4_power))
  
  target <- targetVariance(lm_AR1_power)
  expect_equal(vcovCR(lm_AR1_power, type = "CR4", target = target, inverse_var = TRUE), CR4_power)
  attr(CR4_power, "inverse_var") <- FALSE
  expect_equal(vcovCR(lm_AR1_power, type = "CR4", target = target, inverse_var = FALSE), CR4_power)
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(lm_AR1, vcov = "CR2"))
  expect_true(check_CR(lm_AR1_power, vcov = "CR2"))
  expect_true(check_CR(lm_AR1, vcov = "CR4"))
  expect_true(check_CR(lm_AR1_power, vcov = "CR4"))
})

test_that("getData works.", {
  re_order <- sample(nrow(Ovary))
  egg_scramble <- Ovary[re_order,]
  gls_scramble <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), 
                      data = egg_scramble)
  dat <- getData(gls_scramble)
  expect_identical(egg_scramble, dat)
})


CR_types <- paste0("CR",0:4)

test_that("Order doesn't matter.", {
  re_order <- sample(nrow(Ovary))
  dat_scramble <- Ovary[re_order,]
  lm_scramble <- update(lm_AR1, data = dat_scramble)
  CR_fit <- lapply(CR_types, function(x) vcovCR(lm_AR1, type = x))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(lm_scramble, type = x))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(lm_AR1, vcov = x, test = "All"))
  test_scramble <- lapply(CR_types, function(x) coef_test(lm_scramble, vcov = x, test = "All"))
  expect_equal(test_fit, test_scramble, tolerance = 10^-6)
  
  constraints <- combn(length(coef(lm_AR1)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(lm_AR1, constraints = constraints, vcov = "CR2", test = "All")
  Wald_scramble <- Wald_test(lm_scramble, constraints = constraints, vcov = "CR2", test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})


test_that("clubSandwich works with dropped observations", {
  dat_miss <- Ovary
  dat_miss$follicles[sample.int(nrow(Ovary), size = round(nrow(Ovary) / 10))] <- NA
  lm_dropped <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = dat_miss,
                    correlation = corAR1(form = ~ 1 | Mare), na.action = na.omit)
  lm_complete <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), 
                     data = dat_miss, subset = !is.na(follicles),
                     correlation = corAR1(form = ~ 1 | Mare))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, test = "All"))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, test = "All"))
  expect_identical(test_drop, test_complete)
})


test_that("CR2 is equivalent to Welch t-test for DiD design", {
})
