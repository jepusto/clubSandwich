context("rma.uni objects")
set.seed(20190513)

library(robumeta, quietly=TRUE)
suppressMessages(library(metafor, quietly=TRUE))

data(corrdat)
corr_robu <- robu(effectsize ~ males + college + binge, data = corrdat, 
                   modelweights = "CORR", studynum = studyid,
                   var.eff.size = var)
corrdat$wt <- corr_robu$data.full$r.weights

corr_meta <- rma(effectsize ~ males + college + binge, data = corrdat, 
                 weights = wt, vi = var, method = "FE")

test_that("CR2 t-tests agree with robumeta for correlated effects", {
  robu_CR2 <- vcovCR(corr_meta, cluster = corrdat$studyid, target = 1 / corrdat$wt, type = "CR2")
  expect_true(check_CR(corr_meta, vcov = robu_CR2))
  # expect_true(check_CR(corr_meta, vcov = "CR4", cluster = corrdat$studyid))
  expect_equivalent(as.matrix(robu_CR2), corr_robu$VR.r)
  expect_equivalent(as.matrix(vcovCR(corr_meta, cluster = corrdat$studyid, 
                                     inverse_var = TRUE, type = "CR2")), corr_robu$VR.r)
  
  CR2_ttests <- coef_test(corr_meta, vcov = robu_CR2, test = "Satterthwaite")
  expect_equal(corr_robu$dfs, CR2_ttests$df)
  expect_equal(corr_robu$reg_table$prob, CR2_ttests$p_Satt)
})

data(hierdat)
hier_meta <- rma(effectsize ~ binge + followup + sreport + age, data = hierdat, 
                 vi = var, method = "REML")
hierdat$wt <- with(hier_meta, 1 / (vi + tau2))
hier_robu <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, userweights = wt)

test_that("CR2 t-tests agree with robumeta for user weighting", {
  skip("Skip until robumeta discrepancies resolved.")
  
  robu_CR2_iv <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid)
  robu_CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid,
                         target = hier_robu$data.full$avg.var.eff.size)
  expect_true(check_CR(hier_meta, vcov = robu_CR2_iv))
  # expect_true(check_CR(hier_meta, vcov = "CR4", cluster = hierdat$studyid))
  expect_true(check_CR(hier_meta, vcov = robu_CR2_not))
  # expect_true(check_CR(hier_meta, vcov = "CR4", cluster = hierdat$studyid,
  #                      target = hier_robu$data.full$avg.var.eff.size))

  expect_that(all.equal(hier_robu$VR.r, as.matrix(robu_CR2_iv)), is_a("character"))
  expect_equivalent(hier_robu$VR.r, as.matrix(robu_CR2_not))
  
  CR2_ttests <- coef_test(hier_meta, vcov = robu_CR2_not, test = "Satterthwaite")
  expect_equal(hier_robu$dfs, CR2_ttests$df)
  expect_equal(hier_robu$reg_table$prob, CR2_ttests$p_Satt)
})

test_that("bread works", {
  expect_true(check_bread(corr_meta, cluster = corrdat$studyid, y = corrdat$effectsize))
  X <- model_matrix(corr_meta)
  W <- corr_meta$weights
  V <- corr_meta$vi
  vcov_corr <- crossprod((sqrt(V) * W * X) %*% bread(corr_meta) / nobs(corr_meta))
  attr(vcov_corr, "dimnames") <- attr(vcov(corr_meta), "dimnames")
  expect_equal(vcov(corr_meta), vcov_corr)

  expect_true(check_bread(hier_meta, cluster = hierdat$studyid, y = hierdat$effectsize))
  expect_equal(vcov(hier_meta), bread(hier_meta) / nobs(hier_meta))
})

CR_types <- paste0("CR",0:4)

test_that("order doesn't matter", {
  dat_scramble <- hierdat[sample(nrow(hierdat)),]
  hier_scramble <-  rma(effectsize ~ binge + followup + sreport + age, 
                        data = dat_scramble, vi = var, method = "REML")
  
  CR_fit <- lapply(CR_types, function(x) vcovCR(hier_meta, type = x, cluster = hierdat$studyid))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(hier_scramble, type = x, cluster = dat_scramble$studyid))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(hier_meta, vcov = x, cluster = hierdat$studyid, test = "All", p_values = FALSE))
  test_scramble <- lapply(CR_types, function(x) coef_test(hier_scramble, vcov = x, cluster = dat_scramble$studyid, test = "All", p_values = FALSE))
  expect_equal(test_fit, test_scramble, tolerance = 10^-6)
  
  constraints <- combn(length(coef(hier_scramble)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(hier_meta, constraints = constraints, vcov = "CR2", cluster = hierdat$studyid, test = "All")
  Wald_scramble <- Wald_test(hier_scramble, constraints = constraints, vcov = "CR2", cluster = dat_scramble$studyid, test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})

test_that("clubSandwich works with dropped covariates", {
  dat_miss <- hierdat
  dat_miss$binge[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  dat_miss$followup[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 20))] <- NA
  expect_warning(hier_drop <- rma(effectsize ~ binge + followup + sreport + age, 
                                  data = dat_miss, vi = var, method = "REML"))
  
  
  subset_ind <- with(dat_miss, !is.na(binge) & !is.na(followup))
  hier_complete <- rma(effectsize ~ binge + followup + sreport + age, 
                       subset = !is.na(binge) & !is.na(followup),
                       data = dat_miss, vi = var, method = "REML")
  expect_error(vcovCR(hier_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop_A <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = dat_miss$studyid))
  CR_drop_B <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = hierdat$studyid))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x, cluster = dat_miss$studyid[subset_ind]))
  expect_equal(CR_drop_A, CR_complete)
  expect_equal(CR_drop_B, CR_complete)
  
  test_drop_A <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, cluster = dat_miss$studyid, test = "All", p_values = FALSE))
  test_drop_B <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, cluster = hierdat$studyid, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(hier_complete, vcov = x, cluster = dat_miss$studyid[subset_ind], test = "All", p_values = FALSE))
  expect_equal(test_drop_A, test_complete, tolerance = 10^-6)
  expect_equal(test_drop_B, test_complete, tolerance = 10^-6)
})

test_that("clubSandwich works with missing variances", {
  
  dat_miss <- hierdat
  dat_miss$var[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  expect_warning(hier_drop <- rma(effectsize ~ binge + followup + sreport + age, 
                                  data = dat_miss, vi = var, method = "REML"))
  
  
  subset_ind <- with(dat_miss, !is.na(var))
  hier_complete <- rma(effectsize ~ binge + followup + sreport + age, 
                       subset = !is.na(var),
                       data = dat_miss, vi = var, method = "REML")
  expect_error(vcovCR(hier_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop_A <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = dat_miss$studyid))
  CR_drop_B <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = hierdat$studyid))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x, cluster = dat_miss$studyid[subset_ind]))
  expect_equal(CR_drop_A, CR_complete)
  expect_equal(CR_drop_B, CR_complete)
  
})


test_that("vcovCR options work for CR2", {
  RE_var <- hier_meta$tau2 + hierdat$var
  CR2_iv <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid)
  expect_identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = FALSE)
  attr(CR2_iv, "inverse_var") <- FALSE
  attr(CR2_iv, "target") <- attr(CR2_not, "target")
  expect_equal(CR2_not, CR2_iv)
  expect_identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var), CR2_not)
  expect_identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = hierdat$var), CR2_not))
})

