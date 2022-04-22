context("plm objects - random effects")
set.seed(20190513)

library(nlme, quietly=TRUE)
library(lme4, quietly=TRUE)
library(plm, quietly=TRUE)

data("Grunfeld", package = "plm")
data("Produc", package = "plm")

# grun_re <- plm(inv ~ value + capital, data = Grunfeld, model="random")
# Grunfeld$cluster <- sample(LETTERS[1:10], size = nrow(Grunfeld), replace=TRUE)
# Grunfeld_scramble <- Grunfeld[sample(nrow(Grunfeld)),]

CR_types <- paste0("CR",0:4)

plm_individual <- plm(inv ~ value + capital, data = Grunfeld, model="random")
obj <- plm_individual

test_that("individual effects agree with gls", {
  icc <- with(plm_individual$ercomp, sigma2[["id"]] / (sigma2[["id"]] + sigma2[["idios"]]))
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

  test_plm <- lapply(CR_types, function(x) coef_test(plm_individual, vcov = x, test = "All", p_values = FALSE)[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_individual, vcov = x, test = "All", p_values = FALSE)[,-3])
  expect_equivalent(test_plm, test_gls)
  
})

plm_time <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc, index = c("state","year"), 
                effect = "time", model = "random")

test_that("time effects agree with gls", {
  icc <- with(plm_time$ercomp, sigma2[[2]] / (sigma2[[2]] + sigma2[[1]]))
  gls_time <- gls(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, 
                  correlation = corCompSymm(value = icc, form = ~ 1 | year, fixed = TRUE))
  
  expect_equal(model_matrix(plm_time), model_matrix(gls_time))
  expect_identical(nobs(plm_time), nobs(gls_time))
  expect_equivalent(residuals_CS(plm_time), residuals_CS(gls_time))
  
  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_time, type = x))
  CR_gls <- lapply(CR_types, function(x) vcovCR(gls_time, type = x))
  expect_equivalent(CR_plm, CR_gls)
  
  test_plm <- lapply(CR_types, function(x) coef_test(plm_time, vcov = x, test = "All", p_values = FALSE)[,-3])
  test_gls <- lapply(CR_types, function(x) coef_test(gls_time, vcov = x, test = "All", p_values = FALSE)[,-3])
  expect_equivalent(test_plm, test_gls)
  
})


plm_twoways <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                   data = Produc, index = c("state","year"), 
                   effect = "twoways", model = "random")

test_that("two-way effects throws error", {
  expect_error(vcovCR(plm_twoways, type = "CR2"))
})

plm_nested <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                  data = Produc, index = c("state","year","region"), 
                  effect = "nested", model = "random")

test_that("nested effects agree with lmer", {
  
  Produc_sort_order <- with(Produc, order(region, state))
  
  plm_nested$ercomp$sigma2
  lmer_nested_fit <- lmer(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp + (1 | region / state), 
                          data = Produc)
  theta <- getME(lmer_nested_fit, "theta")
  theta[1:2] <- with(plm_nested$ercomp, sqrt(sigma2[2:3] / sigma2[1]))
  lmer_nested <- update(lmer_nested_fit, start = theta, control = lmerControl(optimizer = NULL))
  
  expect_equivalent(model_matrix(plm_nested), model_matrix(lmer_nested)[Produc_sort_order,])
  expect_identical(nobs(plm_nested), nobs(lmer_nested))
  
  expect_equivalent(targetVariance(plm_nested, cluster = findCluster.plm(plm_nested)),
                    targetVariance(lmer_nested, cluster = Produc$region))
  
  expect_equivalent(weightMatrix(plm_nested, cluster = findCluster.plm(plm_nested)),
                    weightMatrix(lmer_nested, cluster = Produc$region), tol = 1e-6)
  
  CR_plm <- lapply(CR_types, function(x) vcovCR(plm_nested, type = x))
  CR_lmer <- lapply(CR_types, function(x) vcovCR(lmer_nested, type = x))
  expect_equivalent(CR_plm, CR_lmer)
  
  test_plm <- lapply(CR_types, function(x) coef_test(plm_nested, vcov = x, test = "All", p_values = FALSE)[,-3])
  test_lmer <- lapply(CR_types, function(x) coef_test(lmer_nested, vcov = x, test = "All", p_values = FALSE)[,-3])
  compare_ttests(test_plm, test_lmer, tol = 1e-5)
  
})

test_that("bread works", {
  
  expect_true(check_bread(plm_individual, 
                          cluster = findCluster.plm(plm_individual), 
                          y = plm_individual$model$inv))
  expect_true(check_bread(plm_time, 
                          cluster = findCluster.plm(plm_time), 
                          y = plm_time$model$"log(gsp)"))
  expect_true(check_bread(plm_nested, 
                          cluster = findCluster.plm(plm_nested), 
                          y = plm_nested$model$"log(gsp)"))
})

test_that("CR0 and CR1S agree with arellano vcov", {
  
  expect_equivalent(vcovHC(plm_individual, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR0")))
  expect_equivalent(vcovHC(plm_individual, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR1S")))
  
  expect_equivalent(vcovHC(plm_time, method="arellano", type = "HC0", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR0")))
  expect_equivalent(vcovHC(plm_time, method="arellano", type = "sss", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR1S")))

  # Can't replicate vcovHC because plm isn't clustering correctly.  
  # expect_equivalent(vcovHC(plm_nested, method="arellano", type = "HC0", cluster = "group"), 
  #                   as.matrix(vcovCR(plm_nested, type = "CR0")))
  # expect_equivalent(vcovHC(plm_nested, method="arellano", type = "sss", cluster = "group"), 
  #                   as.matrix(vcovCR(plm_nested, type = "CR1S")))
})


test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(plm_individual, type = "CR2")
  expect_equal(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR2"), CR2_iv)
  expect_equal(vcovCR(plm_individual, type = "CR2", inverse_var = TRUE), CR2_iv)
  tgt <- targetVariance(plm_individual, cluster = Grunfeld$firm)
  expect_equivalent(vcovCR(plm_individual, type = "CR2", target = tgt, inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_individual, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_equal(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_equal(vcovCR(plm_individual, type = "CR2", target = tgt), CR2_not)
  expect_equal(vcovCR(plm_individual, type = "CR2", target = tgt, inverse_var = FALSE), CR2_not)
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_individual, type = "CR4")
  expect_equal(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR4"), CR4_iv)
  expect_equal(vcovCR(plm_individual, type = "CR4", inverse_var = TRUE), CR4_iv)
  tgt <- targetVariance(plm_individual, cluster = Grunfeld$firm)
  expect_equivalent(vcovCR(plm_individual, type = "CR4", target = tgt, inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_individual, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_equal(vcovCR(plm_individual, cluster = Grunfeld$firm, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_equal(vcovCR(plm_individual, type = "CR4", target = tgt), CR4_not)
  expect_equal(vcovCR(plm_individual, type = "CR4", target = tgt, inverse_var = FALSE), CR4_not)
})

test_that("CR2 and CR4 are target-unbiased", {
  
  expect_true(check_CR(plm_individual, vcov = "CR2"))
  expect_true(check_CR(plm_individual, vcov = "CR4"))
  
  expect_true(check_CR(plm_time, vcov = "CR2"))
  expect_true(check_CR(plm_time, vcov = "CR4"))
  
  expect_true(check_CR(plm_nested, vcov = "CR2"))
  expect_true(check_CR(plm_nested, vcov = "CR4"))
  
})


test_that("vcovCR works when clustering at a level above the random effects.", {
  
  data("Wages", package = "plm")
  Wages$ID <- rep(1:595, each = 7)
  Wages$period <- rep(1:7, times = 595)
  Wages$Grp <- rep(1:119, each = 7 * 5)
  
  plm_ID <- plm(lwage ~ wks + south + smsa + married + exp, 
                data = Wages, 
                index = c("ID","period","Grp"),
                model="random")
  
  ICC <- with(plm_ID$ercomp, sigma2[2] / sum(sigma2))
  
  gls_ID <- gls(lwage ~ wks + south + smsa + married + exp, 
                data = Wages, 
                correlation = corCompSymm(value = ICC, form = ~ 1 | ID, fixed = TRUE))
  
  plm_vcov_ID <- lapply(CR_types, function(x) vcovCR(plm_ID, type = x))
  gls_vcov_ID <- lapply(CR_types, function(x) vcovCR(gls_ID, type = x))

  expect_equivalent(plm_vcov_ID, gls_vcov_ID)
  
  plm_vcov_grp <- lapply(CR_types, function(x) vcovCR(plm_ID, cluster =  Wages$Grp, type = x))
  plm_vcov_group <- lapply(CR_types, function(x) vcovCR(plm_ID, cluster = "group", type = x))
  expect_equal(plm_vcov_grp, plm_vcov_group)
  
  gls_vcov_grp <- lapply(CR_types, function(x) vcovCR(gls_ID, cluster = Wages$Grp, type = x))
  expect_equivalent(plm_vcov_group, gls_vcov_grp)
  
  plm_ID <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc, index = c("state","year","region"), 
                effect = "individual", model = "random")
  
  ICC <- with(plm_ID$ercomp, sigma2[2] / sum(sigma2))
  
  gls_ID <- gls(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc, 
                correlation = corCompSymm(value = ICC, form = ~ 1 | state, fixed = TRUE))
  
  plm_vcov_grp <- lapply(CR_types, function(x) vcovCR(plm_ID, cluster = Produc$region, type = x))
  plm_vcov_group <- lapply(CR_types, function(x) vcovCR(plm_ID, cluster = "group", type = x))
  expect_equal(plm_vcov_grp, plm_vcov_group)
  
  gls_vcov_grp <- lapply(CR_types, function(x) vcovCR(gls_ID, cluster = Produc$region, type = x))
  expect_equivalent(plm_vcov_grp, gls_vcov_grp)
  
})
