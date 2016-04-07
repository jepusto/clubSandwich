context("rma.mv objects")

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


test_that("withS and withG model specifications agree.", {
  dat_long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  levels(dat_long$group) <- c("exp", "con")
  dat_long$group <- relevel(dat_long$group, ref="con")
  dat_long$esid <- factor(1:nrow(dat_long))
  dat_long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat_long)
  rma_G <- rma.mv(yi, vi, mods = ~ group, random = ~ group | study, struct="CS", data=dat_long)
  rma_S <- rma.mv(yi, vi, mods = ~ group, random = list(~ 1 | esid, ~ 1 | study), data=dat_long)

  CR_G <- lapply(CR_types, function(x) vcovCR(rma_G, type = x))
  CR_S <- lapply(CR_types, function(x) vcovCR(rma_S, type = x))
  expect_equivalent(CR_G, CR_S)
  
  tests_G <- lapply(CR_types, function(x) coef_test(rma_G, vcov = x, test = "All"))
  tests_S <- lapply(CR_types, function(x) coef_test(rma_S, vcov = x, test = "All"))
  expect_equal(tests_G, tests_S, tolerance = 10^-6)
})

test_that("order doesn't matter", {
  dat_scramble <- hierdat[sample(nrow(hierdat)),]
  hier_scramble <-  rma.mv(effectsize ~ binge + followup + sreport + age, 
                           random = list(~ 1 | esid, ~ 1 | studyid),
                           data = dat_scramble, V = var, method = "REML")
  
  CR_fit <- lapply(CR_types, function(x) vcovCR(hier_meta, type = x, cluster = hierdat$studyid))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(hier_scramble, type = x, cluster = dat_scramble$studyid))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(hier_meta, vcov = x, cluster = hierdat$studyid, test = "All"))
  test_scramble <- lapply(CR_types, function(x) coef_test(hier_scramble, vcov = x, cluster = dat_scramble$studyid, test = "All"))
  expect_equal(test_fit, test_scramble, tolerance = 10^-6)
  
  constraints <- combn(length(coef(hier_scramble)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(hier_meta, constraints = constraints, vcov = "CR2", cluster = hierdat$studyid, test = "All")
  Wald_scramble <- Wald_test(hier_scramble, constraints = constraints, vcov = "CR2", cluster = dat_scramble$studyid, test = "All")
  expect_equal(Wald_fit, Wald_scramble)
})

test_that("clubSandwich errors with dropped observations", {
  dat_miss <- hierdat
  dat_miss$binge[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  dat_miss$followup[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 20))] <- NA
  expect_warning(hier_drop <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                                     random = list(~ 1 | esid, ~ 1 | studyid),
                                     data = dat_miss, V = var, method = "REML"))
  expect_error(vcovCR(hier_drop, type = "CR0", cluster = dat_miss$studyid))
  
  hier_complete <- rma.mv(effectsize ~ binge + followup + sreport + age, 
                          random = list(~ 1 | esid, ~ 1 | studyid),
                          subset = !is.na(binge) & !is.na(followup),
                          data = dat_miss, V = var, method = "REML")
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, test = "All"))
  test_complete <- lapply(CR_types, function(x) coef_test(hier_complete, vcov = x, test = "All"))
  expect_equal(test_drop, test_complete, tolerance = 10^-6)
})

test_that("vcovCR options work for CR2", {
  RE_var <- targetVariance(hier_meta, cluster = factor(hierdat$studyid))
  CR2_iv <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid)
  expect_identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, inverse_var = FALSE)
  expect_equal(CR2_not, CR2_iv)
  expect_equivalent(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var), CR2_not)
  expect_equivalent(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = RE_var, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(hier_meta, type = "CR2", cluster = hierdat$studyid, target = hierdat$var), CR2_not))
})
