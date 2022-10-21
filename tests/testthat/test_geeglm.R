context("logit glm objects")
set.seed(202201007)

# library(doBy)
library(geepack)

# For ar1
## With waves
timeorder <- rep(1:5, 6)
tvar      <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each=5)
uuu   <- rep(rnorm(6), each=5)
yvar  <- 1 + 2 * tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)

simdatPerm <- simdat[sample(nrow(simdat)),]
simdatPerm <- simdatPerm[order(simdatPerm$idvar),]
wav <- simdatPerm$timeorder

mod_ar1_wav <- geeglm(yvar ~ tvar, id = idvar, 
                      data = simdatPerm, 
                      corstr = "ar1", waves = wav)

## No waves
data(dietox)
dietox$Cu <- as.factor(dietox$Cu)
mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))
mod_ar1 <- geeglm(mf, data = dietox, 
                  id = Pig, 
                  family = poisson("identity"), 
                  corstr = "ar1")


test_that("bread works", {
  
  expect_true(check_bread(mod_ar1, cluster = dietox$Pig, y = dietox$Weight))
  expect_true(check_bread(mod_ar1_wav, cluster = simdatPerm$idvar, y = simdatPerm$yvar))

})

test_that("vcovCR works for clustering variables higher than id variable.", {
  
})

test_that("vcovCR options work for CR2", {
  
  CR2_iv <- vcovCR(mod_ar1_wav, cluster = mod_ar1_wav$cluster, type = "CR2")
  expect_equal(vcovCR(mod_ar1_wav, cluster = mod_ar1_wav$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(mod_ar1_wav, cluster = mod_ar1_wav$cluster, type = "CR2", 
                          target = targetVariance(mod_ar1_wav, cluster = dat$cluster), 
                          inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(mod_ar1_wav, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)

  CR2_iv <- vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2")
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)
  
  CR2_iv <- vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2")
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                          inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR2_iv)
  attr(CR2_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR2", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR2_iv)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(logit_fit, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(logit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(logit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
  CR4_iv <- vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(sflogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(sflogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
  CR4_iv <- vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4")
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                          inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = TRUE), CR4_iv)
  attr(CR4_iv, "inverse_var") <- FALSE
  expect_equal(vcovCR(plogit_fit, cluster = dat$cluster, type = "CR4", 
                      target = targetVariance(plogit_fit, cluster = dat$cluster), 
                      inverse_var = FALSE), CR4_iv)
  
})

test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(logit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(sflogit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(plogit_fit, vcov = "CR2", cluster = dat$cluster))
  expect_true(check_CR(logit_fit, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(sflogit_fit, vcov = "CR4", cluster = dat$cluster))
  expect_true(check_CR(plogit_fit, vcov = "CR4", cluster = dat$cluster))
})

test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)
  
  HC_types <- paste0("HC", 0:3)
  HC_list <- lapply(HC_types, function(t) vcovHC(logit_fit, type = t))
  
  CR_types <- paste0("CR", 0:3)
  CR_types[2] <- "CR1S"
  CR_list <- lapply(CR_types, function(t) as.matrix(vcovCR(logit_fit, cluster = dat$row, type = t)))
  expect_equal(HC_list, CR_list, tol = 4 * 10^-4)
  
})

CR_types <- paste0("CR", 0:4)

test_that("Order doesn't matter.",{
  check_sort_order(logit_fit, dat = dat, cluster = "cluster")
})


test_that("clubSandwich works with dropped observations", {
  dat_miss <- dat
  dat_miss$X1[sample.int(n, size = round(n / 10))] <- NA
  logit_dropped <- update(logit_fit, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(X1))
  logit_complete <- update(logit_fit, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(logit_dropped, cluster = dat_miss$cluster, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(logit_complete, cluster = dat_complete$cluster, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(logit_dropped, vcov = x, cluster = dat_miss$cluster, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(logit_complete, vcov = x, cluster = dat_complete$cluster, test = "All", p_values = FALSE))
  compare_ttests(test_drop, test_complete)
})


test_that("clubSandwich works with aliased predictors", {
  data(npk, package = "datasets")
  npk_alias <- glm(yield ~ block + N*P*K, data = npk)
  npk_drop <- glm(yield ~ block + N + P + K + N:P + N:K + P:K, data = npk)
  
  CR_alias <- lapply(CR_types[-4], function(x) vcovCR(npk_alias, cluster = npk$block, type = x))
  CR_drop <- lapply(CR_types[-4], function(x) vcovCR(npk_drop, cluster = npk$block, type = x))
  expect_equal(CR_alias, CR_drop)
  
  test_drop <- lapply(CR_types[-4], function(x) coef_test(npk_alias, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite"), coefs = 7:12, p_values = FALSE))
  test_complete <- lapply(CR_types[-4], function(x) coef_test(npk_drop, vcov = x, cluster = npk$block, test = c("z","naive-t","Satterthwaite"), coefs = 7:12, p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


test_that("clubSandwich results are equivalent to geepack", {
  
  skip_if_not_installed("geepack")
  library(geepack)
  
  # check CR0 with logit
  logit_gee <- geeglm(y1 ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, family = "binomial")
  logit_refit <- update(logit_fit, start = coef(logit_gee))
  expect_equal(coef(logit_refit), coef(logit_gee))
  
  V_gee0 <- summary(logit_gee)$cov.scaled
  V_CR0 <- as.matrix(vcovCR(logit_refit, cluster = dat$cluster, type = "CR0"))
  attr(V_gee0, "dimnames") <- attr(V_CR0, "dimnames")
  expect_equal(V_gee0, V_CR0)
  
  # check CR3 with logit
  logit_gee <- geeglm(y1 ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, family = "binomial",
                      std.err = "jack")
  V_gee3 <- summary(logit_gee)$cov.scaled
  V_CR3 <- as.matrix(vcovCR(logit_refit, cluster = dat$cluster, type = "CR3"))
  attr(V_gee3, "dimnames") <- attr(V_CR3, "dimnames")
  expect_equal(V_gee3 * m / (m - 6), V_CR3)
  
  # check CR0 with plogit
  plogit_gee <- geeglm(yp ~ X1 + X2 + X3 + X4, id = cluster, 
                      data = dat, weights = w, family = "binomial")
  plogit_refit <- update(plogit_fit, start = coef(plogit_gee))
  expect_equal(coef(plogit_refit), coef(plogit_gee))
  
  V_gee0 <- summary(plogit_gee)$cov.scaled
  V_CR0 <- as.matrix(vcovCR(plogit_refit, cluster = dat$cluster, type = "CR0"))
  attr(V_gee0, "dimnames") <- attr(V_CR0, "dimnames")
  expect_equal(V_gee0, V_CR0)
  
})
