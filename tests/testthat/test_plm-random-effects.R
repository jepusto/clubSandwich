context("plm objects - random effects")

library(nlme, quietly=TRUE)
library(plm, quietly=TRUE)

data("Grunfeld", package = "plm")
data("Produc", package = "plm")

# grun_re <- plm(inv ~ value + capital, data = Grunfeld, model="random")
# Grunfeld$cluster <- sample(LETTERS[1:10], size = nrow(Grunfeld), replace=TRUE)
# Grunfeld_scramble <- Grunfeld[sample(nrow(Grunfeld)),]

CR_types <- paste0("CR",0:4)

plm_individual <- plm(inv ~ value + capital, data = Grunfeld, model="random")

test_that("individual effects agree with gls", {
  icc <- with(plm_individual$ercomp$sigma2, id / (id + idios))
  gls_individual <- gls(inv ~ value + capital, data = Grunfeld,
                        correlation = corCompSymm(value = icc, form = ~ 1 | firm, fixed=TRUE))
  
  expect_equal(model_matrix(plm_individual), model_matrix(gls_individual))
  expect_identical(nobs(plm_individual), nobs(gls_individual))
  V_ratio <- Map("/", targetVariance(plm_individual, cluster = Grunfeld$firm),
                 targetVariance(gls_individual, cluster = Grunfeld$firm))
  expect_equal(lapply(V_ratio, min), lapply(V_ratio, max))
  expect_equivalent(residuals_CS(plm_individual), residuals_CS(gls_individual))

  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_individual, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_individual, type = x))
  expect_equivalent(CR_plm, CR_gls)

  test_plm <- lapply(CR_types, function(x) coef_test(plm_individual, vcov = x, test = "All")[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_individual, vcov = x, test = "All")[,-3])
  expect_equivalent(test_plm, test_gls)
  
})

plm_time <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc, index = c("state","year"), 
                effect = "time", model = "random")

test_that("time effects agree with gls", {
  icc <- with(plm_time$ercomp$sigma2, id / (id + idios))
  gls_time <- gls(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, 
                  correlation = corCompSymm(value = icc, form = ~ 1 | year, fixed = TRUE))
  
  expect_equal(model_matrix(plm_time), model_matrix(gls_time))
  expect_identical(nobs(plm_time), nobs(gls_time))
  V_ratio <- Map("/", targetVariance(plm_time, cluster = Produc$year),
                 targetVariance(gls_time, cluster = Produc$year))
  expect_equal(lapply(V_ratio, min), lapply(V_ratio, max))
  expect_equivalent(residuals_CS(plm_time), residuals_CS(gls_time))
  
  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_time, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_time, type = x))
  expect_equivalent(CR_plm, CR_gls)
  
  test_plm <- lapply(CR_types, function(x) coef_test(plm_time, vcov = x, test = "All")[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_time, vcov = x, test = "All")[,-3])
  expect_equivalent(test_plm, test_gls)
  
})


plm_twoways <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                   data = Produc, index = c("state","year"), 
                   effect = "twoways", model = "random")

test_that("two-way effects throws error", {
  expect_error(vcovCR(plm_twoways, type = "CR2"))
})

test_that("bread works", {
  expect_true(check_bread(plm_individual, 
                          cluster = findCluster.plm(plm_individual), 
                          y = plm_individual$model$inv))
  expect_equal(vcov(plm_individual), 
               plm_individual$ercomp$sigma2$idios * bread(plm_individual) / v_scale(plm_individual))
  
  expect_true(check_bread(plm_time, 
                          cluster = findCluster.plm(plm_time), 
                          y = plm_time$model$inv))
  expect_equal(vcov(plm_time), 
               plm_time$ercomp$sigma2$idios * bread(plm_time) / v_scale(plm_time))
  
  expect_true(check_bread(plm_twoways, 
                          cluster = findCluster.plm(plm_individual), 
                          y = plm_twoways$model$inv))
  expect_true(check_bread(plm_twoways, 
                          cluster = findCluster.plm(plm_time), 
                          y = plm_twoways$model$inv))
  expect_equal(vcov(plm_twoways), 
               plm_twoways$ercomp$sigma2$idios * bread(plm_twoways) / v_scale(plm_twoways))
})

test_that("CR0 and CR1S agree with arellano vcov", {
  expect_equal(vcovHC(plm_individual, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR0")))
  expect_equal(vcovHC(plm_individual, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR1S")))
  
  expect_equal(vcovHC(plm_time, method="arellano", type = "HC0", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR0")))
  expect_equal(vcovHC(plm_time, method="arellano", type = "sss", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR1S")))
})


test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(plm_individual, type = "CR2")
  expect_identical(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR2"), CR2_iv)
  expect_identical(vcovCR(plm_individual, type = "CR2", inverse_var = TRUE), CR2_iv)
  tgt <- targetVariance(plm_individual, cluster = Grunfeld$firm)
  expect_equivalent(vcovCR(plm_individual, type = "CR2", target = tgt, inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_individual, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_identical(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_identical(vcovCR(plm_individual, type = "CR2", target = tgt), CR2_not)
  expect_identical(vcovCR(plm_individual, type = "CR2", target = tgt, inverse_var = FALSE), CR2_not)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_individual, type = "CR4")
  expect_identical(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR4"), CR4_iv)
  expect_identical(vcovCR(plm_individual, type = "CR4", inverse_var = TRUE), CR4_iv)
  tgt <- targetVariance(plm_individual, cluster = Grunfeld$firm)
  expect_equivalent(vcovCR(plm_individual, type = "CR4", target = tgt, inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_individual, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_identical(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_identical(vcovCR(plm_individual, type = "CR4", target = tgt), CR4_not)
  expect_identical(vcovCR(plm_individual, type = "CR4", target = tgt, inverse_var = FALSE), CR4_not)
})

test_that("CR2 and CR4 are target-unbiased", {
  
  expect_true(check_CR(plm_individual, vcov = "CR2"))
  expect_true(check_CR(plm_individual, vcov = "CR4"))
  
  expect_true(check_CR(plm_time, vcov = "CR2"))
  expect_true(check_CR(plm_time, vcov = "CR4"))

})
