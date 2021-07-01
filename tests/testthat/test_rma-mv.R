context("rma.mv objects")
set.seed(20190513)

library(robumeta, quietly=TRUE)
suppressMessages(library(metafor, quietly=TRUE))

data(corrdat)
corr_robu <- robu(effectsize ~ males + college + binge, data = corrdat, 
                   modelweights = "CORR", studynum = studyid,
                   var.eff.size = var)
corrdat$wt <- corr_robu$data.full$r.weights
corr_meta <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                    V = var, W = wt, method = "FE")

test_that("CR2 t-tests agree with robumeta for correlated effects", {
  
  robu_CR2 <- vcovCR(corr_meta, cluster = corrdat$studyid, target = 1 / corrdat$wt, type = "CR2")
  expect_true(check_CR(corr_meta, vcov = robu_CR2))
  # expect_true(check_CR(corr_meta, vcov = "CR4", cluster = corrdat$studyid))
  expect_equivalent(as.matrix(robu_CR2), corr_robu$VR.r)
  expect_that(all.equal(as.matrix(vcovCR(corr_meta, cluster = corrdat$studyid, 
                                     inverse_var = TRUE, type = "CR2")), corr_robu$VR.r),
                    is_a("character"))
  
  CR2_ttests <- coef_test(corr_meta, vcov = robu_CR2, test = "Satterthwaite")
  expect_equal(corr_robu$dfs, CR2_ttests$df)
  expect_equal(corr_robu$reg_table$prob, CR2_ttests$p_Satt)
})

data(hierdat)
hier_meta <- rma.mv(effectsize ~ binge + followup + sreport + age, data = hierdat, 
                    random = list(~ 1 | esid, ~ 1 | studyid),
                    V = var, method = "REML")
hier_robu <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, modelweights = "HIER")

test_that("CR2 t-tests do not exactly agree with robumeta for hierarchical weighting", {
  
  robu_CR2_iv <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid)
  robu_CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid,
                         target = hier_robu$data.full$avg.var.eff.size)
  expect_true(check_CR(hier_meta, vcov = robu_CR2_iv))
  # expect_true(check_CR(hier_meta, vcov = "CR4"))
  expect_true(check_CR(hier_meta, vcov = robu_CR2_not))
  # expect_true(check_CR(hier_meta, vcov = "CR4", 
  #                      target = hier_robu$data.full$avg.var.eff.size))
  
  expect_that(all.equal(hier_robu$VR.r, as.matrix(robu_CR2_iv), check.attributes=FALSE), is_a("character"))
  expect_that(all.equal(hier_robu$VR.r, as.matrix(robu_CR2_not), check.attributes=FALSE), is_a("character"))
  
  CR2_ttests <- coef_test(hier_meta, vcov = robu_CR2_not, test = "Satterthwaite")
  expect_that(all.equal(hier_robu$dfs, CR2_ttests$df), is_a("character"))
  expect_that(all.equal(hier_robu$reg_table$prob, CR2_ttests$p_Satt), is_a("character"))
})

CR_types <- paste0("CR",0:4)

dat_long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
levels(dat_long$group) <- c("exp", "con")
dat_long$group <- relevel(dat_long$group, ref="con")
dat_long$esid <- factor(1:nrow(dat_long))
dat_long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat_long)
rma_G <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="CS", data=dat_long)
rma_S <- rma.mv(yi, vi, mods = ~ group, random = list(~ 1 | esid, ~ 1 | study), data=dat_long)

test_that("withS and withG model specifications agree.", {
  CR_G <- lapply(CR_types, function(x) vcovCR(rma_G, type = x))
  CR_S <- lapply(CR_types, function(x) vcovCR(rma_S, type = x))
  expect_equivalent(CR_G, CR_S)
  
  tests_G <- lapply(CR_types, function(x) coef_test(rma_G, vcov = x, test = "All", p_values = FALSE))
  tests_S <- lapply(CR_types, function(x) coef_test(rma_S, vcov = x, test = "All", p_values = FALSE))
  expect_equal(tests_G, tests_S, tolerance = 10^-6)
})

test_that("bread works", {
  expect_true(check_bread(corr_meta, cluster = corrdat$studyid, y = corrdat$effectsize))
  X <- model_matrix(corr_meta)
  W <- corr_meta$W
  V <- corr_meta$vi
  vcov_corr <- bread(corr_meta) %*% t(X) %*% W %*% (V * W) %*% X %*% bread(corr_meta) / nobs(corr_meta)^2
  attr(vcov_corr, "dimnames") <- attr(vcov(corr_meta), "dimnames")
  expect_equal(vcov(corr_meta), vcov_corr)
  
  expect_true(check_bread(hier_meta, cluster = hierdat$studyid, y = hierdat$effectsize))
  expect_equal(vcov(hier_meta), bread(hier_meta) / nobs(hier_meta))
  
  expect_true(check_bread(rma_G, cluster = dat_long$study, y = dat_long$yi))
  expect_equal(vcov(rma_G), bread(rma_G) / nobs(rma_G))
  
  expect_true(check_bread(rma_S, cluster = dat_long$study, y = dat_long$yi))
  expect_equal(vcov(rma_S), bread(rma_S) / nobs(rma_S))
})

test_that("order doesn't matter", {

  skip_on_cran()
  
  check_sort_order(hier_meta, hierdat)

})

test_that("clubSandwich works with dropped covariates", {
  dat_miss <- hierdat
  dat_miss$binge[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  dat_miss$followup[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 20))] <- NA
  expect_warning(hier_drop <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                                     random = list(~ 1 | esid, ~ 1 | studyid),
                                     data = dat_miss, V = var, method = "REML"))
  
  hier_complete <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                          random = list(~ 1 | esid, ~ 1 | studyid),
                          subset = !is.na(binge) & !is.na(followup),
                          data = dat_miss, V = var, method = "REML")
  
  expect_error(vcovCR(hier_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop_A <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x))
  CR_drop_B <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = dat_miss$studyid))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x))
  expect_equal(CR_drop_A, CR_complete)
  expect_equal(CR_drop_B, CR_complete)
  
  test_drop_A <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, test = "All", p_values = FALSE))
  test_drop_B <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, cluster = dat_miss$studyid, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(hier_complete, vcov = x, test = "All", p_values = FALSE))
  compare_ttests(test_drop_A, test_complete)
  compare_ttests(test_drop_B, test_complete)
  
})

test_that("clubSandwich works with missing diagonal variances", {
  
  dat_miss <- hierdat
  dat_miss$var[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  expect_warning(hier_drop <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                                     random = list(~ 1 | esid, ~ 1 | studyid),
                                     data = dat_miss, V = var, method = "REML"))
  
  hier_complete <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                          random = list(~ 1 | esid, ~ 1 | studyid),
                          subset = !is.na(var),
                          data = dat_miss, V = var, method = "REML")
  
  expect_error(vcovCR(hier_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop_A <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x))
  CR_drop_B <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x, cluster = dat_miss$studyid))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x))
  expect_equal(CR_drop_A, CR_complete)
  expect_equal(CR_drop_B, CR_complete)
  
})

test_that("clubSandwich works with missing vcov matrix", {
  
  skip_if(packageVersion("metafor") < 2.1)
  
  dat_miss <- corrdat
  dat_miss$var[sample.int(nrow(corrdat), size = round(nrow(corrdat) / 10))] <- NA
  
  V_missing <- impute_covariance_matrix(dat_miss$var, cluster = dat_miss$studyid, r = 0.8)
  
  expect_warning(corr_drop <- rma.mv(effectsize ~ males + college + binge, 
                                     random = ~ 1 | studyid, 
                                     V = V_missing, data = dat_miss))
  
  corr_complete <- rma.mv(effectsize ~ males + college + binge,
                          random = ~ 1 | studyid, 
                          subset = !is.na(var),
                          data = dat_miss, V = V_missing)
  
  expect_error(vcovCR(corr_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(corr_drop, cluster = dat_miss$studyid, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(corr_complete, type = x))
  expect_equal(CR_drop, CR_complete, tol = 10^-5)
  
  
  # V_complete <- impute_covariance_matrix(corrdat$var, cluster = corrdat$studyid, r = 0.8)
  # W_missing <- lapply(V_complete, function(x) chol2inv(chol(x)))
  # 
  # corr_drop <- rma.mv(effectsize ~ males + college + binge, 
  #                     random = ~ 1 | studyid, 
  #                     V = V_complete, W = bldiag(W_missing), 
  #                     data = dat_miss)
  # 
  # corr_complete <- rma.mv(effectsize ~ males + college + binge,
  #                         random = ~ 1 | studyid, 
  #                         V = V_complete, W = bldiag(W_missing),
  #                         data = dat_miss, subset = !is.na(var))
  # 
  # expect_error(vcovCR(corr_complete, type = "CR0", cluster = dat_miss$studyid))
  # 
  # CR_drop <- lapply(CR_types, function(x) vcovCR(corr_drop, type = x))
  # CR_complete <- lapply(CR_types, function(x) vcovCR(corr_complete, type = x))
  # expect_equal(CR_drop, CR_complete)
  
})


test_that("vcovCR options work for CR2", {
  RE_var <- targetVariance(hier_meta, cluster = factor(hierdat$studyid))
  CR2_iv <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid)
  expect_equal(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = FALSE)
  expect_equal(CR2_not, CR2_iv)
  expect_equivalent(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var), CR2_not)
  expect_equivalent(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = hierdat$var), CR2_not))
})

test_that("clubSandwich works with complicated random effects specifications.", {
  
  skip_on_cran()
  
  data(oswald2013, package = "robumeta")
  
  oswald2013 <- within(oswald2013, {
    V = (1 - R^2)^2 / (N - 3)
    SSID = paste(Study, "sample",Sample.ID)
    ESID = 1:nrow(oswald2013)
  })
  
  SS_lab <- unique(oswald2013$SSID)
  n_SS <- length(SS_lab)
  R_mat <- 0.4 + 0.6 * diag(nrow = n_SS)
  colnames(R_mat) <- rownames(R_mat) <- SS_lab
  
  m1 <- rma.mv(
    R ~ 0 + IAT.Focus + Crit.Cat, V = V,
    data = oswald2013,
    random = list(~ 1 | Study, ~ 1 | SSID, ~ 1 | ESID)
  )
  
  m2 <- update(m1, 
               random = list(~ IAT.Focus | Study, ~ 1 | SSID, ~ 1 | ESID),
               struct = c("UN"))
  
  m3 <- update(m1, 
               random = list(~ 1 | Study, ~ IAT.Focus | SSID, ~ 1 | ESID),
               struct = c("UN","UN"))
  
  m4 <- update(m1, 
               random = list(~ 1 | Study, ~ 1 | SSID, ~ IAT.Focus | ESID),
               struct = c("DIAG"))
  
  m5 <- update(m1, 
               random = list(~ IAT.Focus | Study, ~ IAT.Focus | SSID),
               struct = c("UN","UN"))
  
  m6 <- update(m1, 
               random = list(~ IAT.Focus | Study, ~ IAT.Focus | SSID, ~ 1 | ESID),
               struct = c("UN","UN"))
  
  m7 <- update(m5, struct = c("CS","CS"))
  m8 <- update(m5, struct = c("HCS","HCS"))
  m9 <- update(m5, struct = c("UN","CS"))
  m10 <- update(m5, struct = c("CS","UN"))
  
  m11 <- rma.mv(
    R ~ 0 + IAT.Focus + Crit.ID, V = V,
    data = oswald2013,
    random = list(~ 1  + Crit.ID | Study),
    struct = c("GEN")
  )
  
  m12 <- update(m11, random = list(~ 1  + Crit.ID | Study, ~ 1 | SSID))
  m13 <- update(m11, random = list(~ 1  + Crit.ID | Study, ~ IAT.Focus | SSID),
                struct = c("GEN","UN"))
  m14 <- update(m11, random = list(~ IAT.Focus | Study, ~ 1 + Crit.ID | SSID),
                struct = c("UN","GEN"))
  
  mod_list <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
  os_cluster <- factor(oswald2013$Study)
  
  obj <- m6
  struct <- parse_structure(obj)
  findCluster.rma.mv(obj)
  
  cluster_list <- lapply(mod_list, findCluster.rma.mv)
  lapply(cluster_list, expect_equal, expected = os_cluster)
  
  bread_checks <- sapply(mod_list, check_bread, cluster = oswald2013$Study, y = oswald2013$R)
  expect_true(all(bread_checks))

  CR_checks <- sapply(mod_list, check_CR, vcov = "CR2")
  expect_true(all(CR_checks))

  m11 <- update(m1, R = list(SSID = R_mat))
  expect_error(findCluster.rma.mv(m11))
  
})

test_that("clubSandwich works for correlated hierarchical effects model.", {
  
  skip_on_cran()
  
  V_mat <- impute_covariance_matrix(vi = corrdat$var, 
                                    cluster = corrdat$studyid,
                                    r = 0.7,
                                    smooth_vi = TRUE)
  
  CHE_es <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                   V = V_mat, random = ~ 1 | esid)
  CHE_study <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                      V = V_mat, random = ~ 1 | studyid)
  CHE_studyes <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                        V = V_mat, random = ~ 1 | studyid / esid)
  CHE_esstudy <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                        V = V_mat, random = ~ 1 | esid/ studyid)
  CHE_study_es <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                         V = V_mat, random = list(~ 1 | studyid,  ~ 1 | esid))
  CHE_es_study <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                         V = V_mat, random = list(~ 1 | esid, ~ 1 | studyid))
  
  mods <- list(es = CHE_es, study = CHE_study, 
               studyes = CHE_studyes, esstudy = CHE_esstudy,
               study_es = CHE_study_es, es_study = CHE_es_study)
  
  clusters <- lapply(mods, findCluster.rma.mv)
  
  expect_equal(clusters$study, clusters$studyes)
  expect_equal(clusters$es, clusters$esstudy)
  expect_equal(clusters$studyes, clusters$study_es)
  expect_equal(clusters$study_es, clusters$es_study)
  
  V_CR2s <- lapply(mods, vcovCR, type = "CR2")
  V_CR2s_clust <- mapply(vcovCR, mods, clusters, type = "CR2", SIMPLIFY = FALSE)
  
  expect_equal(V_CR2s, V_CR2s_clust)
  expect_equal(V_CR2s$studyes, V_CR2s$study_es)
  expect_equal(V_CR2s$studyes, V_CR2s$es_study)
  
})