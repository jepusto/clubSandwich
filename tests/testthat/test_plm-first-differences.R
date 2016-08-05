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

# obj <- plm_FD
# vcov <- vcovCR(obj, type = "CR4")
# type <- "CR4"
# dim(model_matrix(obj))
# 
# index <- attr(model.frame(obj),"index")
# cluster <- switch(obj$args$effect,
#                   individual = index[[1]],
#                   time = index[[2]])
# sort_order <- get_index_order(obj)
# cluster <- cluster[sort_order]
# 
# if (obj$args$model=="fd") {
#   cluster <- cluster[index[[2]] != levels(index[[2]])[1]]
# }
# 
# target <- NULL
# inverse_var <- is.null(target)
# obj$na.action <- attr(obj$model, "na.action")
# 
# vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)

test_that("CR0 and CR1S agree with arellano vcov", {
  expect_equal(vcovHC(plm_FD, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_FD, type = "CR0")))
  expect_equal(vcovHC(plm_FD, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_FD, type = "CR1S")))
  
  expect_equal(vcovHC(plm_FD, method="arellano", type = "HC0", cluster = "time"),
               as.matrix(vcovCR(plm_FD, cluster = "time", type = "CR0")))
  expect_equal(vcovHC(plm_FD, method="arellano", type = "sss", cluster = "time"),
               as.matrix(vcovCR(plm_FD, cluster = "time", type = "CR1S")))
  
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
  expect_equal(CR_individual, HC_individual)
  
})
