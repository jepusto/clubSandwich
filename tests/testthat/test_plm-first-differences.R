context("plm objects - first differences model")

library(plm, quietly=TRUE)

data(Fatalities, package = "AER")
Fatalities <- within(Fatalities, {
  frate <- 10000 * fatal / pop
  drinkagec <- cut(drinkage, breaks = 18:22, include.lowest = TRUE, right = FALSE)
  drinkagec <- relevel(drinkagec, ref = 4)
})

plm_FD <- plm(frate ~ beertax + drinkagec + miles + unemp + log(income),
              data = Fatalities, index = c("state", "year"), 
              model = "fd")

n_obs <- nobs(plm_FD)
target <- with(Fatalities, 1 / pop[year != levels(year)[1]])

test_that("bread works", {
  y <- na.omit(diff(plm_FD$model$frate))
  cluster <- findCluster.plm(plm_FD)
  expect_true(check_bread(plm_FD, cluster = cluster, y = y))
  sigma_sq <- with(plm_FD, sum(residuals^2) / df.residual) 
  expect_equal(vcov(plm_FD), bread(plm_FD) * sigma_sq / v_scale(plm_FD))
})

test_that("CR0 and CR1S agree with arellano vcov", {
  
  expect_equal(vcovHC(plm_FD, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_FD, type = "CR0")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_FD, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_FD, type = "CR1S")), check.attributes = FALSE)
  
  X <- model_matrix(plm_FD)
  e <- residuals(plm_FD)
  index <- attr(model.frame(plm_FD), "index")
  cluster <- index[[2]]
  cluster <- cluster[index[[2]] != levels(index[[2]])[1]]
  estmats <- sapply(split.data.frame(e * X, cluster, drop = TRUE), colSums)
  meat <- tcrossprod(estmats)
  bread <- chol2inv(chol(crossprod(X)))
  vcov_time <- bread %*% meat %*% bread
  attr(vcov_time, "dimnames") <- attr(meat, "dimnames")
  
  expect_equal(vcov_time, as.matrix(vcovCR(plm_FD, cluster = "time", type = "CR0")))
  
  baloney <- tcrossprod(estmats[,-6])
  vcov_baloney <- bread %*% baloney %*% bread
  attr(vcov_baloney, "dimnames") <- attr(baloney, "dimnames")
  expect_equal(vcov_baloney, 
               vcovHC(plm_FD, method="arellano", type = "HC0", cluster = "time"), check.attributes = FALSE)
})

test_that("vcovCR options work for CR2", {
  
  CR2_iv <- vcovCR(plm_FD, type = "CR2")
  expect_identical(vcovCR(plm_FD, cluster = Fatalities$state, type = "CR2"), CR2_iv)
  expect_identical(vcovCR(plm_FD, type = "CR2", inverse_var = TRUE), CR2_iv)
  
  expect_identical(vcovCR(plm_FD, type = "CR2", 
                          target = rep(1, n_obs), 
                          inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_FD, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_identical(vcovCR(plm_FD, cluster = Fatalities$state, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_identical(vcovCR(plm_FD, type = "CR2", target = rep(1, n_obs)), CR2_not)
  expect_identical(vcovCR(plm_FD, type = "CR2", target = rep(1, n_obs), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(plm_FD, type = "CR2", target = target), CR2_not))
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_FD, type = "CR4")
  expect_identical(vcovCR(plm_FD, cluster = Fatalities$state, type = "CR4"), CR4_iv)
  expect_identical(vcovCR(plm_FD, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_identical(vcovCR(plm_FD, type = "CR4", target = rep(1, n_obs), inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_FD, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_identical(vcovCR(plm_FD, cluster = Fatalities$state, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_identical(vcovCR(plm_FD, type = "CR4", target = rep(1, n_obs)), CR4_not)
  expect_identical(vcovCR(plm_FD, type = "CR4", target = rep(1, n_obs), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(plm_FD, type = "CR4", target = target), CR4_not))
})

test_that("CR2 is target-unbiased", {
  
  expect_true(check_CR(plm_FD, vcov = "CR2"))
  expect_true(check_CR(plm_FD, vcov = "CR2", inverse_var = FALSE))
  
  expect_true(check_CR(plm_FD, cluster = "time", vcov = "CR2"))
  expect_true(check_CR(plm_FD, cluster = "time", vcov = "CR2", inverse_var = FALSE))

})


test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  CR_types <- paste0("CR",c(0,2))
  HC_types <- paste0("HC",c(0,2))

  CR_individual <- lapply(CR_types, function(t) 
    as.matrix(vcovCR(plm_FD, cluster = 1:nrow(Fatalities), type = t)))
  HC_individual <- lapply(HC_types, function(t) 
    vcovHC(plm_FD, method = "white1", type = t))
  expect_equal(CR_individual, HC_individual, check.attributes = FALSE)
  
})

