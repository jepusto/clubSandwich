context("plm objects")

library(nlme)
library(plm)

data("Grunfeld", package = "plm")
data("Produc", package = "plm")

# grun_re <- plm(inv ~ value + capital, data = Grunfeld, model="random")
# Grunfeld$cluster <- sample(LETTERS[1:10], size = nrow(Grunfeld), replace=TRUE)
# Grunfeld_scramble <- Grunfeld[sample(nrow(Grunfeld)),]

CR_types <- paste0("CR",0:4)

test_that("individual effects agree with gls", {
  plm_individual <- plm(inv ~ value + capital, data = Grunfeld, model="random")

  icc <- with(plm_individual$ercomp$sigma2, id / (id + idios))
  gls_individual <- gls(inv ~ value + capital, data = Grunfeld,
                        correlation = corCompSymm(value = icc, form = ~ 1 | firm))

  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_individual, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_individual, type = x))
  expect_less_than(max(mapply(function(x, y) max(abs(range(x / y) - 1)), 
                              x = CR_plm, y = CR_gls)), .06)

  test_plm <- lapply(CR_types, function(x) coef_test(plm_individual, vcov = x, test = "All")[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_individual, vcov = x, test = "All")[,-3])
  expect_less_than(max(mapply(function(x, y) max(abs(range(x / y) - 1)), 
                              x = test_plm, y = test_gls)), .10)
  
})

test_that("time effects agree with gls", {
  plm_time <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, index = c("state","year"), 
                  effect = "time", model = "random")
  icc <- with(plm_time$ercomp$sigma2, id / (id + idios))
  gls_time <- gls(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, 
                  correlation = corCompSymm(value = icc, form = ~ 1 | year))
  
  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_time, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_time, type = x))
  expect_less_than(max(mapply(function(x, y) max(abs(range(x / y) - 1)), 
                              x = CR_plm, y = CR_gls)), .06)
  
  test_plm <- lapply(CR_types, function(x) coef_test(plm_time, vcov = x, test = "All")[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_time, vcov = x, test = "All")[,-3])
  expect_less_than(max(mapply(function(x, y) max(abs(range(x / y) - 1)), 
                              x = test_plm, y = test_gls)), .10)
})

test_that("two-way effects throws error", {
  plm_twoways <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                          data = Produc, index = c("state","year"), 
                          effect = "twoways", model = "random")
  
  expect_error(vcovCR(plm_twoways, type = "CR2"))
})

# test cluster specification
# test target matrix specification
# test inverse_var detection