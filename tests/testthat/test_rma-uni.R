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
  
  skip_on_cran()
  
  check_sort_order(hier_meta, hierdat, cluster = "studyid")
  
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
  compare_ttests(test_drop_A, test_complete)
  compare_ttests(test_drop_B, test_complete)
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
  expect_equal(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = FALSE)
  attr(CR2_iv, "inverse_var") <- FALSE
  attr(CR2_iv, "target") <- attr(CR2_not, "target")
  expect_equal(CR2_not, CR2_iv)
  expect_equal(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var), CR2_not)
  expect_equal(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = hierdat$var), CR2_not))
})


test_that("vcovCR works with intercept-only model and user-specified weights.", {
  
  dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  dat$wt <- sample(1:3, size = nrow(dat), replace = TRUE)
  
  res <- rma(yi, vi, weights = wt, data=dat)
  
  meta_rob <- robust(res, cluster=dat$trial)
  club_rob <- coef_test(res, vcov="CR1", cluster=dat$trial, test = "naive-t")
  expect_equal(meta_rob$se, club_rob$SE)
  expect_equal(meta_rob$zval, club_rob$tstat)
  expect_equal(meta_rob$dfs, club_rob$df_t)
  expect_equal(meta_rob$pval, club_rob$p_t)
  
  expect_true(check_CR(res, vcov = "CR2", cluster = dat$trial))
  test_uni <- coef_test(res, vcov="CR2", cluster=dat$trial, test = "All")
  
  
  res <- rma.mv(yi, vi, W = wt, random = ~ 1 | trial, data=dat) 
  meta_rob <- robust(res, cluster=dat$trial)
  club_rob <- coef_test(res, vcov="CR1", test = "naive-t")
  expect_equal(meta_rob$se, club_rob$SE)
  expect_equal(meta_rob$zval, club_rob$tstat)
  expect_equal(meta_rob$dfs, club_rob$df_t)
  expect_equal(meta_rob$pval, club_rob$p_t)
  
  expect_true(check_CR(res, vcov = "CR2"))
  test_mv <- coef_test(res, vcov="CR2", test = "All")
  
  expect_equal(test_uni, test_mv, tolerance = 10^-5)

  V_club <- vcovCR(res, type = "CR2")
  k <- res$k
  yi <- res$yi
  wi <- diag(res$W)
  W <- sum(wi)
  wi <- wi / W
  vi <- diag(res$M)
  V <- sum(vi)
  ei <- residuals_CS(res)
  M <- sum(wi^2 * vi)
  ai <- 1 / sqrt(1 - 2 * wi + M / vi)
  V_hand <- sum(wi^2 * ai^2 * ei^2)
  expect_equal(V_hand, as.numeric(V_club))
  
  pi_theta_pj <- diag(vi) - tcrossprod(rep(1,k), wi * vi) - tcrossprod(wi * vi, rep(1, k)) + M
  df <- M^2 / sum(tcrossprod(ai^2 * wi^2) * (pi_theta_pj^2))
  
  expect_equal(Inf, test_uni$df_z)
  expect_equal(k - 1, test_uni$df_t)
  expect_equal(df, test_uni$df_Satt, tolerance = 10^-5)
  
})
