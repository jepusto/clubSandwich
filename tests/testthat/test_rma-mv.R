context("rma.mv objects")
set.seed(20190513)

skip_if_not_installed("robumeta")
skip_if_not_installed("metafor")

CR_types <- paste0("CR",0:4)

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
  expect_equal(tests_G, tests_S, tolerance = 1e-6)
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
  
  skip_if(packageVersion("metafor") < "2.1")
  
  dat_miss <- corrdat
  dat_miss$var[sample.int(nrow(corrdat), size = round(nrow(corrdat) / 10))] <- NA
  
  V_missing <- impute_covariance_matrix(dat_miss$var, cluster = dat_miss$studyid, r = 0.8)
  
  expect_warning(corr_drop <- rma.mv(effectsize ~ males + college + binge, 
                                     random = ~ 1 | studyid, 
                                     V = V_missing, data = dat_miss,
                                     sparse = TRUE))
  
  corr_complete <- rma.mv(effectsize ~ males + college + binge,
                          random = ~ 1 | studyid, 
                          subset = !is.na(var),
                          data = dat_miss, V = V_missing,
                          sparse = TRUE)
  
  expect_error(vcovCR(corr_complete, type = "CR0", cluster = dat_miss$studyid))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(corr_drop, cluster = dat_miss$studyid, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(corr_complete, type = x))
  expect_equal(CR_drop, CR_complete, tolerance = 1e-5)
  
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
    random = list(~ 1 | Study, ~ 1 | SSID, ~ 1 | ESID),
    sparse = TRUE
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
    struct = c("GEN"),
    sparse = TRUE
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

test_that("clubSandwich works for random slopes model.", {
  
  # example from https://wviechtb.github.io/metadat/reference/dat.obrien2003.html
  dat <- dat.obrien2003
  dat$bmicent <- dat$bmi - ave(dat$bmi, dat$study)
  dat <- escalc(measure="PR", xi=cases, ni=total, data=dat)
  dat$yi <- dat$yi*100
  dat$vi <- dat$vi*100^2
  res <- rma.mv(yi, vi, mods = ~ bmicent, 
                random = ~ bmicent | study, struct="GEN", 
                data=dat,
                sparse = TRUE)
  
  cl <- findCluster.rma.mv(res)
  
  expect_true(check_bread(res, cluster = cl, y = dat$yi))
  expect_true(check_CR(res, vcov = "CR2"))
  
})

test_that("clubSandwich works for correlated hierarchical effects model.", {
  
  skip_on_cran()
  
  V_mat <- impute_covariance_matrix(vi = corrdat$var, 
                                    cluster = corrdat$studyid,
                                    r = 0.7,
                                    smooth_vi = TRUE)
  
  CHE_es <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                   V = V_mat, random = ~ 1 | esid,
                   sparse = TRUE)
  CHE_study <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                      V = V_mat, random = ~ 1 | studyid,
                      sparse = TRUE)
  CHE_studyes <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                        V = V_mat, random = ~ 1 | studyid / esid,
                        sparse = TRUE)
  CHE_esstudy <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                        V = V_mat, random = ~ 1 | esid/ studyid,
                        sparse = TRUE)
  CHE_study_es <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                         V = V_mat, random = list(~ 1 | studyid,  ~ 1 | esid),
                         sparse = TRUE)
  CHE_es_study <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                         V = V_mat, random = list(~ 1 | esid, ~ 1 | studyid),
                         sparse = TRUE)
  
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

test_that("vcovCR errors when there is only one cluster.", {
  
  dat <- data.frame(
    study = "study1", # study number
    est = runif(5, 0.1, 0.6), # R-squared values
    se = runif(5, 0.005, 0.025), # standard errors of R-squared values
    es_id = 1:5 # effect size ID
  )
  
  v_mat <- impute_covariance_matrix(dat$se^2, cluster = dat$study, r = 0.8)
  
  # working model in metafor
  expect_warning(
    res <- rma.mv(yi = est, V = v_mat, random = ~ 1 | study / es_id, data = dat)
  )

  single_cluster_error_msg <- "Cluster-robust variance estimation will not work when the data only includes a single cluster."

  expect_error(
    vcovCR(res, type = "CR0"), single_cluster_error_msg
  )

  expect_error(
    conf_int(res, vcov = "CR1"), single_cluster_error_msg
  )
  
  expect_error(
    coef_test(res, vcov = "CR2"), single_cluster_error_msg
  )
  
  expect_error(
    Wald_test(res, constraints = constrain_zero(1), vcov = "CR3"),
    single_cluster_error_msg
  )
  
  expect_error(
    vcovCR(res, cluster = dat$es_id),
    "Random effects are not nested within clustering variable."
  )
  
})


test_that("clubSandwich works when random effects variable has missing levels.",{
  
  dat <- dat.konstantopoulos2011
  dat$district_fac <- factor(dat$district)
  dat$district_fac_plus <- factor(dat$district, levels = c(levels(dat$district_fac), 1000, 10000))
  
  mlma_fac <- rma.mv(yi ~ year, V = vi, 
                     random = ~ 1 | district_fac / study,
                     data = dat,
                     sparse = TRUE)
  
  implicit_fac <- coef_test(mlma_fac, vcov = "CR2")
  explicit_fac <- coef_test(mlma_fac, vcov = "CR2", cluster = dat$district_fac)
  expect_equal(implicit_fac, explicit_fac)

  mlma_plus <- rma.mv(yi ~ year, V = vi, 
                     random = ~ 1 | district_fac_plus / study,
                     data = dat,
                     sparse = TRUE)
  
  implicit_plus <- coef_test(mlma_plus, vcov = "CR2")
  explicit_plus <- coef_test(mlma_plus, vcov = "CR2", cluster = dat$district_fac_plus)
  expect_equal(implicit_plus, explicit_plus)
  expect_equal(implicit_fac, implicit_plus)
  expect_equal(implicit_fac, explicit_plus)
  
  mlma_num <- rma.mv(yi ~ year, V = vi, 
                     random = ~ 1 | district / study,
                     data = dat,
                     sparse = TRUE)
  
  implicit_num <- coef_test(mlma_num, vcov = "CR2")
  explicit_num <- coef_test(mlma_num, vcov = "CR2", cluster = dat$district)
  expect_equal(implicit_num, explicit_num)
  expect_equal(implicit_fac, implicit_num)
  expect_equal(implicit_fac, explicit_num)
  
})

Vmat <- with(corrdat, impute_covariance_matrix(vi = var, cluster = studyid, r = 0.8))
corr_meta <- rma.mv(effectsize ~ males + college + binge, data = corrdat, 
                    V = Vmat, random = ~ 1 | studyid)

test_that("clubSandwich agrees with metafor::robust() for CR0.", {
  
  test_CR0 <- coef_test(corr_meta, vcov = "CR0", test = "All")
  meta_CR0 <- robust(corr_meta, cluster = corrdat$studyid, adjust = FALSE)
  rob_CR0 <- coef_test(meta_CR0, vcov = "CR0", test = "All")
  expect_equal(test_CR0$SE, meta_CR0$se)
  expect_equal(test_CR0$df_tp, rep(meta_CR0$df, length(test_CR0$df_tp)))
  expect_equal(test_CR0$p_tp, meta_CR0$pval, tolerance = 1e-5)
  compare_ttests(rob_CR0, test_CR0, tol = 1e-6)

  club_F_CR0 <- Wald_test(corr_meta, constraints = constrain_zero(2:4), 
                          vcov = "CR0", test = "Naive-Fp")
  rob_F_CR0 <- Wald_test(meta_CR0, constraints = constrain_zero(2:4), 
                         vcov = "CR0", test = "Naive-Fp")
  expect_equal(club_F_CR0$Fstat, meta_CR0$QM)
  expect_equal(club_F_CR0$df_num, meta_CR0$QMdf[1])
  expect_equal(club_F_CR0$df_denom, meta_CR0$QMdf[2])
  expect_equal(club_F_CR0$p_val, meta_CR0$QMp, tolerance = 1e-5)
  compare_Waldtests(club_F_CR0, rob_F_CR0, tol = 1e-5)

})

test_that("clubSandwich agrees with metafor::robust() for CR1p.", {
  
  test_CR1 <- coef_test(corr_meta, vcov = "CR1p", test = "All")
  meta_CR1 <- robust(corr_meta, cluster = corrdat$studyid, adjust = TRUE)
  rob_CR1 <- coef_test(meta_CR1, vcov = "CR1p", test = "All")
  expect_equal(test_CR1$SE, meta_CR1$se)
  expect_equal(test_CR1$df_tp, rep(meta_CR1$df, length(test_CR1$df_tp)))
  expect_equal(test_CR1$p_tp, meta_CR1$pval, tolerance = 1e-5)
  compare_ttests(rob_CR1, test_CR1, tol = 1e-5)
  
  club_F_CR1 <- Wald_test(corr_meta, constraints = constrain_zero(2:4), 
                          vcov = "CR1p", test = "Naive-Fp")
  rob_F_CR1 <- Wald_test(meta_CR1, constraints = constrain_zero(2:4), 
                         vcov = "CR1p", test = "Naive-Fp")
  expect_equal(club_F_CR1$Fstat, meta_CR1$QM)
  expect_equal(club_F_CR1$df_num, meta_CR1$QMdf[1])
  expect_equal(club_F_CR1$df_denom, meta_CR1$QMdf[2])
  expect_equal(club_F_CR1$p_val, meta_CR1$QMp, tolerance = 1e-5)
  compare_Waldtests(club_F_CR1, rob_F_CR1, tol = 1e-5)
  
})

test_that("clubSandwich agrees with metafor::robust() for CR2.", {
  
  skip_if(packageVersion('metafor') < "3.1.31")
  test_CR2 <- coef_test(corr_meta, vcov = "CR2", test = "All")
  meta_CR2 <- robust(corr_meta, cluster = corrdat$studyid, clubSandwich = TRUE)
  rob_CR2 <- coef_test(meta_CR2, vcov = "CR2", test = "All")
  expect_equal(test_CR2$SE, meta_CR2$se)
  compare_ttests(rob_CR2, test_CR2, tol = 1e-5)
  
  club_F_CR2 <- Wald_test(corr_meta, constraints = constrain_zero(2:4), 
                          vcov = "CR2", test = "All")
  rob_F_CR2 <- Wald_test(meta_CR2, constraints = constrain_zero(2:4), 
                         vcov = "CR2", test = "All")
  expect_equal(subset(club_F_CR2, test == "HTZ")$Fstat, meta_CR2$QM)
  expect_equal(subset(club_F_CR2, test == "HTZ")$df_num, meta_CR2$QMdf[1])
  expect_equal(subset(club_F_CR2, test == "HTZ")$df_denom, meta_CR2$QMdf[2])
  expect_equal(subset(club_F_CR2, test == "HTZ")$p_val, meta_CR2$QMp, tolerance = 1e-5)
  compare_Waldtests(club_F_CR2, rob_F_CR2, tol = 1e-5)
})

test_that("clubSandwich methods work on robust.rma objects.", {
  
  hier_robust <- robust(hier_meta, cluster = hierdat$studyid, adjust = TRUE)
  
  expect_equal(residuals_CS(hier_meta), residuals_CS(hier_robust))
  expect_equal(coef_CS(hier_meta), coef_CS(hier_robust))
  expect_equal(model_matrix(hier_meta), model_matrix(hier_robust))
  expect_equal(bread(hier_meta), bread(hier_robust))
  expect_equal(v_scale(hier_meta), v_scale(hier_robust))
  expect_equal(targetVariance(hier_meta, cluster = hierdat$studyid), 
               targetVariance(hier_robust, cluster = hierdat$studyid))
  expect_equal(weightMatrix(hier_meta, cluster = hierdat$studyid), 
               weightMatrix(hier_robust, cluster = hierdat$studyid))
  
  hier_club <- robust(hier_meta, cluster = hierdat$studyid, adjust = FALSE, clubSandwich = TRUE)
  
  expect_equal(residuals_CS(hier_meta), residuals_CS(hier_club))
  expect_equal(coef_CS(hier_meta), coef_CS(hier_club))
  expect_equal(model_matrix(hier_meta), model_matrix(hier_club))
  expect_equal(bread(hier_meta), bread(hier_club))
  expect_equal(v_scale(hier_meta), v_scale(hier_club))
  expect_equal(targetVariance(hier_meta, cluster = hierdat$studyid), 
               targetVariance(hier_club, cluster = hierdat$studyid))
  expect_equal(weightMatrix(hier_meta, cluster = hierdat$studyid), 
               weightMatrix(hier_club, cluster = hierdat$studyid))
  
})

test_that("clubSandwich works with user-weighted rma.mv objects.", {
  
  data("oswald2013", package = "robumeta")
  oswald2013$yi <- atanh(oswald2013$R)
  oswald2013$vi <- 1 / (oswald2013$N - 3)
  oswald2013$esID <- 1:nrow(oswald2013)
  oswald2013$wt <- 1 + rpois(nrow(oswald2013), lambda = 1)
  
  V <- impute_covariance_matrix(vi = oswald2013$vi, cluster = oswald2013$Study, r = 0.4)
  
  mod_wt1 <- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                    V = V, W = wt,
                    random = ~ 1 | Study,
                    data = oswald2013,
                    sparse = TRUE)
  
  W_mat <- impute_covariance_matrix(vi = oswald2013$wt, cluster = oswald2013$Study, r = 0, return_list = FALSE)
  
  mod_wt2 <- rma.mv(yi ~ 0 + Crit.Cat + Crit.Domain + IAT.Focus + Scoring,
                    V = V, W = W_mat,
                    random = ~ 1 | Study,
                    data = oswald2013,
                    sparse = TRUE)
  
  
  vcovs_1 <- lapply(CR_types, function(x) vcovCR(mod_wt1, type = x))
  vcovs_2 <- lapply(CR_types, function(x) vcovCR(mod_wt2, type = x))

  coef_test_wt1 <- lapply(CR_types, function(x) 
    coef_test(mod_wt1, vcov = x, test = "All")
  )
  
  coef_test_wt2 <- lapply(CR_types, function(x) 
    coef_test(mod_wt2, vcov = x, test = "All")
  )
  
  Wald_test_wt1 <- lapply(CR_types, function(x) 
    Wald_test(mod_wt1, 
              constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
              vcov = x, 
              test = "All")
  )
  
  Wald_test_wt2 <- lapply(CR_types, function(x) 
    Wald_test(mod_wt2, 
              constraints = constrain_equal("Crit.Cat", reg_ex = TRUE),
              vcov = x, 
              test = "All")
  )
  
  expect_equal(vcovs_1, vcovs_2, tolerance = 1e-5)
  compare_ttests(coef_test_wt1, coef_test_wt2, tol = 1e-5)
  compare_Waldtests(Wald_test_wt1, Wald_test_wt2, tol = 1e-5)
  
  for (i in seq_along(vcovs_1)) {
    expect_s3_class(vcovs_1[[i]], "vcovCR")
    expect_s3_class(vcovs_2[[i]], "vcovCR")
    expect_s3_class(coef_test_wt1[[i]], "coef_test_clubSandwich")
    expect_s3_class(coef_test_wt2[[i]], "coef_test_clubSandwich")
    expect_s3_class(Wald_test_wt1[[i]], "Wald_test_clubSandwich")
    expect_s3_class(Wald_test_wt2[[i]], "Wald_test_clubSandwich")
  }
  
  
})
