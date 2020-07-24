context("ivreg objects")
set.seed(20190513)

library(zoo, quietly=TRUE)
library(AER, quietly=TRUE)

data("CigarettesSW", package = "AER")

Cigs <- within(CigarettesSW, {
  rprice <- price/cpi
  rincome <- income/population/cpi
  tdiff <- (taxs - tax)/cpi
})

CR_types <- paste0("CR",0:4)

obj_un <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
             data = Cigs)
obj_wt <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
            data = Cigs, 
            weights = population)

X <- model.matrix(obj_wt, component = "regressors")
Z <- model.matrix(obj_wt, component = "instruments")
y <- log(CigarettesSW$packs)
w <- weights(obj_wt)

test_that("Basic calculations from ivreg agree for unweighted model.", {
  XZ <- model.matrix(obj_un, component = "projected")
  ZtZ_inv <- chol2inv(chol(t(Z) %*% Z))
  XZ_check <- Z %*% ZtZ_inv %*% t(Z) %*% X
  
  expect_equal(XZ, XZ_check, check.attributes=FALSE)
  expect_equal(coef(obj_un), lm.fit(XZ, y)$coefficients)
  expect_equal (bread(obj_un), chol2inv(chol(t(XZ) %*% XZ)) * nobs(obj_un), check.attributes=FALSE)
  
  hii <- diag(X %*% chol2inv(chol(t(XZ) %*% XZ)) %*% t(XZ))
  expect_equal(hatvalues(obj_un), hii)
  
  r <- as.vector(y - X %*% coef(obj_un))
  expect_equal(r, as.vector(residuals_CS(obj_un)))
})

test_that("Basic calculations from ivreg agree for weighted model.", {
  XZ <- model.matrix(obj_wt, component = "projected")
  ZwZ_inv <- chol2inv(chol(t(Z) %*% (w * Z)))
  XZ_check <- Z %*% ZwZ_inv %*% t(Z) %*% (w * X)
  
  expect_equal(XZ, XZ_check, check.attributes=FALSE)
  expect_equal(coef(obj_wt), lm.wfit(XZ, y, w)$coefficients)
  expect_equal(bread(obj_wt), chol2inv(chol(t(XZ) %*% (w * XZ))) * nobs(obj_wt), check.attributes=FALSE)
  
  hii <- diag(X%*% chol2inv(chol(t(XZ) %*% (w * XZ))) %*% t(w * XZ))
  expect_false(all(hatvalues(obj_wt) == hii)) # does not agree because hatvalues doesn't work with weighting
  
  r <- as.vector(y - X %*% coef(obj_wt))
  expect_equal(r, as.vector(residuals_CS(obj_wt)))
})

test_that("bread works", {
  
  expect_true(check_bread(obj_un, cluster = Cigs$state, y = log(Cigs$packs)))
  tsls_vcov <- bread(obj_un) * summary(obj_un)$sigma^2 / v_scale(obj_un)
  expect_equal(vcov(obj_un), tsls_vcov)
  
  expect_true(check_bread(obj_wt, cluster = Cigs$state, y = log(Cigs$packs)))
  wtsls_vcov <- bread(obj_wt) * summary(obj_wt)$sigma^2 / v_scale(obj_wt)
  expect_equal(vcov(obj_wt), wtsls_vcov)
})


test_that("vcovCR options don't matter for CR0", {
  expect_error(vcovCR(obj_un, type = "CR0"))
  CR0 <- vcovCR(obj_un, cluster = Cigs$state, type = "CR0")
  expect_output(print(CR0))
  attr(CR0, "target") <- NULL
  attr(CR0, "inverse_var") <- NULL
  CR0_A <- vcovCR(obj_un, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population)
  attr(CR0_A, "target") <- NULL
  attr(CR0_A, "inverse_var") <- NULL
  expect_identical(CR0_A, CR0)
  CR0_B <- vcovCR(obj_un, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population, inverse_var = FALSE)
  attr(CR0_B, "target") <- NULL
  attr(CR0_B, "inverse_var") <- NULL
  expect_identical(CR0_A, CR0)
  
  expect_error(vcovCR(obj_un, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population, inverse_var = TRUE))

  wCR0 <- vcovCR(obj_wt, cluster = Cigs$state, type = "CR0")
  attr(wCR0, "target") <- NULL
  attr(wCR0, "inverse_var") <- NULL
  wCR0_A <- vcovCR(obj_wt, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population)
  attr(wCR0_A, "target") <- NULL
  attr(wCR0_A, "inverse_var") <- NULL
  expect_identical(wCR0_A, wCR0)
  wCR0_B <- vcovCR(obj_wt, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population, inverse_var = FALSE)
  attr(wCR0_B, "target") <- NULL
  attr(wCR0_B, "inverse_var") <- NULL
  expect_identical(wCR0_B, wCR0)
  
  expect_error(vcovCR(obj_wt, cluster = Cigs$state, type = "CR0", target = 1 / Cigs$population, inverse_var = TRUE))
})

test_that("vcovCR options work for CR2", {
  
  CR2_iv <- vcovCR(obj_un, cluster = Cigs$state, type = "CR2")
  expect_equal(vcovCR(obj_un, cluster = Cigs$state, type = "CR2", target = rep(1, nobs(obj_un))), CR2_iv)
  expect_false(identical(vcovCR(obj_un, cluster = Cigs$state, type = "CR2", target = 1 / Cigs$population), CR2_iv))
  
  wCR2_id <- vcovCR(obj_wt, cluster = Cigs$state, type = "CR2")
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR2", inverse_var = FALSE), wCR2_id)
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR2", target = rep(1, nobs(obj_un))), wCR2_id)
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR2", target = rep(1, nobs(obj_un)), inverse_var = FALSE), wCR2_id)
  
})

test_that("vcovCR options work for CR4", {
  
  CR4_not <- vcovCR(obj_un, cluster = Cigs$state, type = "CR4")
  expect_identical(vcovCR(obj_un, cluster = Cigs$state, type = "CR4", target = rep(1, nobs(obj_un))), CR4_not)
  expect_identical(vcovCR(obj_un, cluster = Cigs$state, type = "CR4", target = rep(1, nobs(obj_un)), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(obj_un, cluster = Cigs$state, type = "CR4", target = 1 / Cigs$population), CR4_not))
  
  wCR4_id <- vcovCR(obj_wt, cluster = Cigs$state, type = "CR4")
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR4", inverse_var = FALSE), wCR4_id)
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR4", target = rep(1, nobs(obj_wt))), wCR4_id)
  expect_identical(vcovCR(obj_wt, cluster = Cigs$state, type = "CR4", target = rep(1, nobs(obj_wt)), inverse_var = FALSE), wCR4_id)
  
})


test_that("CR2 is target-unbiased", {
  expect_true(check_CR(obj_un, vcov = "CR2", cluster = Cigs$state))
  expect_true(check_CR(obj_wt, vcov = "CR2", cluster = Cigs$state))
})

test_that("CR4 is target-unbiased", {
  skip("Need to understand target-unbiasedness for ivreg objects.")
  expect_true(check_CR(obj_un, vcov = "CR4", cluster = Cigs$state))
  expect_true(check_CR(obj_wt, vcov = "CR4", cluster = Cigs$state))
})

test_that("vcovCR is equivalent to vcovHC (with HC0 or HC1) when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)
  CR0 <- vcovCR(obj_un, cluster = 1:nobs(obj_un), type = "CR0")
  expect_equal(vcovHC(obj_un, type = "HC0"), as.matrix(CR0))
  CR1 <- vcovCR(obj_un, cluster = 1:nobs(obj_un), type = "CR1S")
  expect_equal(vcovHC(obj_un, type = "HC1"), as.matrix(CR1))
  CR2 <- vcovCR(obj_un, cluster = 1:nobs(obj_un), type = "CR2")
  expect_false(all(vcovHC(obj_un, type = "HC2") == as.matrix(CR2)))
  CR3 <- vcovCR(obj_un, cluster = 1:nobs(obj_un), type = "CR3")
  expect_false(all(vcovHC(obj_un, type = "HC3") == as.matrix(CR3)))
})

test_that("Order doesn't matter.",{
  
  check_sort_order(obj_wt, Cigs, "state")
  
})

test_that("clubSandwich works with dropped observations", {
  dat_miss <- Cigs
  dat_miss$rincome[sample.int(nrow(Cigs), size = round(nrow(Cigs) / 10))] <- NA
  iv_dropped <- update(obj_un, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(rincome))
  iv_complete <- update(obj_un, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(iv_dropped, cluster = dat_miss$state, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(iv_complete, cluster = dat_complete$state, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(iv_dropped, vcov = x, cluster = dat_miss$state, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(iv_complete, vcov = x, cluster = dat_complete$state, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})



test_that("weight scale doesn't matter", {
  
  iv_fit_w <- update(obj_un, weights = rep(4, nobs(obj_un)))
  
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(obj_un, cluster = Cigs$state, type = x))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(iv_fit_w, cluster = Cigs$state, type = x))
  
  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix), 
               tol = 5 * 10^-7)  
  
  target <- 1 + rpois(nrow(Cigs), lambda = 8)
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(obj_un, cluster = Cigs$state, type = x, target = target))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(iv_fit_w, cluster = Cigs$state, type = x, target = target * 15))
  
  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix), 
               tol = 5 * 10^-7)  
  
})


