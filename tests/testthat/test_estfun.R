context("estfun objects")

library(zoo, quietly=TRUE)
library(AER, quietly=TRUE)

CR_types <- paste0("CR",0:4)

data("CigarettesSW", package = "AER")

Cigs <- within(CigarettesSW, {
  rprice <- price/cpi
  rincome <- income/population/cpi
  tdiff <- (taxs - tax)/cpi
})

obj_un <- ivreg(log(packs) ~ log(rprice) + log(rincome) + I(tax/cpi) | log(rincome) + tdiff + I(tax/cpi),
                data = Cigs)
obj_wt <- ivreg(log(packs) ~ log(rprice) + log(rincome) + I(tax/cpi) | log(rincome) + tdiff + I(tax/cpi),
                data = Cigs, 
                weights = population)

red_form_un <- lm(log(packs) ~ log(rincome) + I(tax/cpi) + tdiff, data = Cigs)
red_form_wt <- lm(log(packs) ~ log(rincome) + I(tax/cpi) + tdiff, data = Cigs, weights = population)
stage1_un <- lm(log(rprice) ~ log(rincome) + I(tax/cpi) + tdiff, data = Cigs)
stage1_wt <- lm(log(rprice) ~ log(rincome) + I(tax/cpi) + tdiff, data = Cigs, weights = population)

test_that("estfun works for lm.", {
  
  V_CR <- lapply(CR_types, function(type) as.matrix(vcovCR(obj = red_form_un, cluster = Cigs$state, type = type)))
  e_CR <- lapply(CR_types, function(type) vcovCR(obj = red_form_un, cluster = Cigs$state, type = type, form = "estfun"))
  expect_equal(lapply(e_CR, tcrossprod), V_CR)
  
  V_CR <- lapply(CR_types, function(type) as.matrix(vcovCR(obj = red_form_wt, cluster = Cigs$state, type = type)))
  e_CR <- lapply(CR_types, function(type) vcovCR(obj = red_form_wt, cluster = Cigs$state, type = type, form = "estfun"))
  expect_equal(lapply(e_CR, tcrossprod), V_CR)

  V_CR <- lapply(CR_types, function(type) as.matrix(vcovCR(obj = stage1_un, cluster = Cigs$state, type = type)))
  e_CR <- lapply(CR_types, function(type) vcovCR(obj = stage1_un, cluster = Cigs$state, type = type, form = "estfun"))
  expect_equal(lapply(e_CR, tcrossprod), V_CR)

  V_CR <- lapply(CR_types, function(type) as.matrix(vcovCR(obj = stage1_wt, cluster = Cigs$state, type = type)))
  e_CR <- lapply(CR_types, function(type) vcovCR(obj = stage1_wt, cluster = Cigs$state, type = type, form = "estfun"))
  expect_equal(lapply(e_CR, tcrossprod), V_CR)
  
})

test_that("stacked estimating equations are equivalent to 2SLS.", {
  
  e_CR <- lapply(CR_types, function(type) vcovCR(obj = red_form_un, cluster = Cigs$state, type = type, form = "estfun"))
  f_CR <- lapply(CR_types, function(type) vcovCR(obj = stage1_un, cluster = Cigs$state, type = type, form = "estfun"))
  
  V_CR_stack <- mapply(function(x, y) tcrossprod(rbind(x, y)), x = e_CR, y = f_CR, SIMPLIFY = FALSE)
  gamma <- coef(stage1_un)["tdiff"]
  beta <- coef(red_form_un)["tdiff"]
  delta <- beta / gamma
  
  V_CR_2SLS <- lapply(CR_types, function(type) vcovCR(obj = obj_un, cluster = Cigs$state, type = type))
  V_CR_2SLS <- sapply(V_CR_2SLS, function(x) diag(x)["log(rprice)"])
  V_delta <- sapply(V_CR_stack, function(x) sum(x[c(4,8), c(4,8)] * tcrossprod(c(1,-delta))) / gamma^2)
})