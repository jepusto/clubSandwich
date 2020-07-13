context("robu objects")
set.seed(20190513)

library(robumeta, quietly=TRUE)
data(corrdat)

corr_large <- robu(effectsize ~ males + college + binge, data = corrdat, 
                   modelweights = "CORR", studynum = studyid,
                   var.eff.size = var, small = FALSE)

test_that("CR0 z-tests agree with robumeta for correlated effects", {
  p <- length(coef_CS(corr_large))
  N <- corr_large$N
  robu_CR0 <- vcovCR(corr_large, type = "CR0")
  ztests <- coef_test(corr_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(corr_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(corr_large$reg_table$SE, ztests$SE)
  expect_equal(with(corr_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
                    ztests$p_z)
})

corr_small <- robu(effectsize ~ males + college + binge, data = corrdat, 
                   modelweights = "CORR", studynum = studyid,
                   var.eff.size = var)

test_that("CR2 t-tests agree with robumeta for correlated effects", {
  robu_CR2 <- vcovCR(corr_small, type = "CR2")
  expect_true(check_CR(corr_small, vcov = robu_CR2))
  # expect_true(check_CR(corr_small, vcov = "CR4"))

  CR2_ttests <- coef_test(corr_small, vcov = robu_CR2, test = "Satterthwaite")
  expect_equivalent(corr_small$VR.r, as.matrix(robu_CR2))
  expect_equal(corr_small$dfs, CR2_ttests$df)
  expect_equal(corr_small$reg_table$prob, CR2_ttests$p_Satt)
})

data(hierdat)

hier_large <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, modelweights = "HIER", small = FALSE)

test_that("CR0 z-tests agree with robumeta for hierarchical effects", {
  p <- length(coef_CS(hier_large))
  N <- hier_large$N
  robu_CR0 <- vcovCR(hier_large, type = "CR0")
  ztests <- coef_test(hier_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(hier_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(hier_large$reg_table$SE, ztests$SE)
  expect_equal(with(hier_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
               ztests$p_z)
})

hier_small <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, modelweights = "HIER")

test_that("CR2 t-tests agree with robumeta for hierarchical effects", {
  robu_CR2 <- vcovCR(hier_small, type = "CR2")
  expect_true(check_CR(hier_small, vcov = robu_CR2))
  # expect_true(check_CR(hier_small, vcov = "CR4"))

  CR2_ttests <- coef_test(hier_small, vcov = robu_CR2, test = "Satterthwaite")
  
  expect_equivalent(hier_small$VR.r, as.matrix(robu_CR2))
  expect_equal(hier_small$dfs, CR2_ttests$df)
  expect_equal(hier_small$reg_table$prob, CR2_ttests$p_Satt)
})

hierdat$user_wt <- 1 + rpois(nrow(hierdat), lambda = 3)

user_large <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, userweights = user_wt, small = FALSE)

test_that("CR0 z-tests agree with robumeta for user weighting", {
  p <- length(coef_CS(user_large))
  N <- user_large$N
  robu_CR0 <- vcovCR(user_large, type = "CR0")
  ztests <- coef_test(user_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(user_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(user_large$reg_table$SE, ztests$SE)
  expect_equal(with(user_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
               ztests$p_z)
})

user_small <- robu(effectsize ~ binge + followup + sreport + age,
                   data = hierdat, studynum = studyid,
                   var.eff.size = var, userweights = user_wt)

test_that("CR2 t-tests agree with robumeta for user weighting", {
  
  user_lm <- lm(effectsize ~ binge + followup + sreport + age, data = hierdat,
                weights = user_wt)
  expect_equivalent(coef_CS(user_lm), coef(user_lm))
  
  robu_CR2 <- vcovCR(user_small, type = "CR2")
  expect_true(check_CR(user_small, vcov = robu_CR2))
  # expect_true(check_CR(user_small, vcov = "CR4"))
  
  expect_equivalent(user_small$VR.r, as.matrix(robu_CR2))
  
  target <- user_small$data.full$avg.var.eff.size
  
  lm_CR2 <- vcovCR(user_lm, cluster = hierdat$studyid, type = "CR2", target = target)
  expect_equivalent(robu_CR2, lm_CR2)
  
  CR2_ttests <- coef_test(user_small, vcov = robu_CR2, test = "Satterthwaite", p_values = FALSE)
#   expect_equal(user_small$dfs, CR2_ttests$df)
#   expect_equal(user_small$reg_table$prob, CR2_ttests$p_Satt)
  lm_CR2_ttests <- coef_test(user_lm, vcov = "CR2", 
                             cluster = hierdat$studyid, 
                             target = user_small$data.full$avg.var.eff.size,
                             test = "Satterthwaite", p_values = FALSE)
  compare_ttests(CR2_ttests, lm_CR2_ttests)
})


test_that("bread works", {
  vcov_corr_large <- with(corr_large, chol2inv(chol(crossprod(Xreg, data.full$r.weights * Xreg))))
  expect_equal(vcov_corr_large, bread(corr_large) / v_scale(corr_large))
  
  vcov_corr_small <- with(corr_small, chol2inv(chol(crossprod(Xreg, data.full$r.weights * Xreg))))
  expect_equal(vcov_corr_small, bread(corr_small) / v_scale(corr_small))
  
  vcov_hier_large <- with(hier_large, chol2inv(chol(crossprod(Xreg, data.full$r.weights * Xreg))))
  expect_equal(vcov_hier_large, bread(hier_large) / v_scale(hier_large))
  
  vcov_hier_small <- with(hier_small, chol2inv(chol(crossprod(Xreg, data.full$r.weights * Xreg))))
  expect_equal(vcov_hier_small, bread(hier_small) / v_scale(hier_small))
  
  vcov_user_large <- with(user_large, chol2inv(chol(crossprod(Xreg, data.full$userweights * Xreg))))
  expect_equal(vcov_user_large, bread(user_large) / v_scale(user_large))
  
  vcov_user_small <- with(user_small, chol2inv(chol(crossprod(Xreg, data.full$userweights * Xreg))))
  expect_equal(vcov_user_small, bread(user_small) / v_scale(user_small))
})


data(dropoutPrevention)

test_that("dropoutPrevention tests replicate Tipton & Pustejovsky (2015) - full sample", {
  skip_on_cran()
  m3_hier <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                  + outcome + evaluator_independence
                  + male_pct + white_pct + average_age
                  + implementation_quality + program_site + duration + service_hrs, 
                  data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")
  
  m3_hier_CR2 <- vcovCR(m3_hier, cluster = dropoutPrevention$studyID, type = "CR2")
  expect_true(check_CR(m3_hier, vcov = m3_hier_CR2))
  # expect_true(check_CR(m3_hier, vcov = "CR4"))
  CR2_ttests <- coef_test(m3_hier, vcov = m3_hier_CR2, test = "Satterthwaite")
  
  expect_equivalent(m3_hier$VR.r, as.matrix(m3_hier_CR2))
  expect_equal(m3_hier$dfs, CR2_ttests$df)
  expect_equal(m3_hier$reg_table$prob, CR2_ttests$p_Satt)

  contrast_list <- list("Study design" = 2:3, 
                        "Outcome measure" = 7:9,
                        "Evaluator independence" = 10:12,
                        "Implmentation quality" = 16:17,
                        "Program format" = 18:20)
  
  dropout_tests <- Wald_test(m3_hier, constraints = constrain_zero(contrast_list), 
                             vcov = m3_hier_CR2, test = c("Naive-F","HTZ"))
  
  Fstat_club <- sapply(dropout_tests, function(x) x$Fstat)
  attr(Fstat_club, "dimnames") <- NULL
  Fstat_paper <- matrix(c(0.23, 0.22, 0.91, 0.84, 3.11, 2.78, 14.15, 13.78, 3.85, 3.65), nrow = 2)
  expect_equivalent(Fstat_paper, round(Fstat_club, 2))
  
  df_club <- sapply(dropout_tests, function(x) x$df_denom[2])
  df_paper <- c(42.9, 21.5, 16.8, 36.9, 37.5)
  attr(df_club, "names") <- NULL
  expect_equivalent(df_paper, round(df_club, 1))
})

test_that("dropoutPrevention tests replicate Tipton & Pustejovsky (2015) - reduced sample", {
  skip_on_cran()
  dp_subset <- subset(dropoutPrevention, big_study==TRUE)
  m3_hier <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                  + outcome + evaluator_independence
                  + male_pct + white_pct + average_age
                  + implementation_quality + program_site + duration + service_hrs, 
                  data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")

  
  m3_hier_CR2 <- vcovCR(m3_hier, cluster = dp_subset$studyID, type = "CR2")
  expect_true(check_CR(m3_hier, vcov = m3_hier_CR2))
  
  CR2_ttests <- coef_test(m3_hier, vcov = m3_hier_CR2, test = "Satterthwaite")
  
  expect_equivalent(m3_hier$VR.r, as.matrix(m3_hier_CR2))
  expect_equal(m3_hier$dfs, CR2_ttests$df)
  expect_equal(m3_hier$reg_table$prob, CR2_ttests$p_Satt)
  
  contrast_list <- list("Study design" = 2:3, 
                        "Outcome measure" = 7:9,
                        "Evaluator independence" = 10:11,
                        "Implmentation quality" = 15:16,
                        "Program format" = 17:19)
  
  dropout_tests <- Wald_test(m3_hier, constraints = constrain_zero(contrast_list), 
                             vcov = "CR2", test = c("Naive-F","HTZ"))
  
  Fstat_club <- sapply(dropout_tests, function(x) x$Fstat)
  Fstat_paper <- matrix(c(3.19, 2.93, 1.05, 0.84, 0.32, 0.26, 4.02, 3.69, 1.19, 0.98), nrow = 2)
  attr(Fstat_club, "dimnames") <- NULL
  expect_equivalent(Fstat_paper, round(Fstat_club, 2))
  
  df_club <- sapply(dropout_tests, function(x) x$df_denom[2])
  df_paper <- c(11.0, 7.7, 4.6, 11.0, 9.1)
  attr(df_club, "names") <- NULL
  expect_equivalent(df_paper, round(df_club, 1))
})

CR_types <- paste0("CR",0:4)

test_that("order doesn't matter", {
  
  skip_on_cran()
  
  dat_scramble <- corrdat[sample(nrow(corrdat)),]
  corr_scramble <-  robu(effectsize ~ males + college + binge, data = dat_scramble, 
                         modelweights = "CORR", studynum = studyid,
                         var.eff.size = var)
  
  CR_fit <- lapply(CR_types, function(x) vcovCR(corr_small, type = x))
  CR_scramble <- lapply(CR_types, function(x) vcovCR(corr_scramble, type = x))
  expect_equivalent(CR_fit, CR_scramble)
  
  test_fit <- lapply(CR_types, function(x) coef_test(corr_small, vcov = x, test = "All", p_values = FALSE))
  test_scramble <- lapply(CR_types, function(x) coef_test(corr_scramble, vcov = x, test = "All", p_values = FALSE))
  compare_ttests(test_fit, test_scramble)
  
  constraints <- combn(length(coef_CS(corr_small)), 2, simplify = FALSE)
  Wald_fit <- Wald_test(corr_small, constraints = constrain_zero(constraints), vcov = "CR2", test = "All")
  Wald_scramble <- Wald_test(corr_scramble, constraints = constrain_zero(constraints), vcov = "CR2", test = "All")
  compare_Waldtests(Wald_fit, Wald_scramble)
})

test_that("clubSandwich works with dropped observations", {
  dat_miss <- hierdat
  dat_miss$binge[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 10))] <- NA
  dat_miss$followup[sample.int(nrow(hierdat), size = round(nrow(hierdat) / 20))] <- NA
  hier_drop <- robu(effectsize ~ binge + followup + sreport + age,
                     data = dat_miss, studynum = studyid,
                     var.eff.size = var, modelweights = "HIER")
  
  dat_complete <- subset(dat_miss, !is.na(binge) & !is.na(followup))
  hier_complete <- robu(effectsize ~ binge + followup + sreport + age,
                    data = dat_complete, studynum = studyid,
                    var.eff.size = var, modelweights = "HIER")
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(hier_drop, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(hier_complete, type = x))
  expect_identical(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(hier_drop, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(hier_complete, vcov = x, test = "All", p_values = FALSE))
  expect_identical(test_drop, test_complete)
})

test_that("vcovCR options work for CR2", {
  dp_subset <- subset(dropoutPrevention, big_study==TRUE)
  m3_hier <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                  + outcome + evaluator_independence
                  + male_pct + white_pct + average_age
                  + implementation_quality + program_site + duration + service_hrs, 
                  data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")
  
  iv <- mean(m3_hier$data.full$r.weights) / m3_hier$data.full$r.weights
  CR2_iv <- vcovCR(m3_hier, type = "CR2")
  expect_identical(vcovCR(m3_hier, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(m3_hier, type = "CR2", target = iv, inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(m3_hier, type = "CR2", inverse_var = FALSE)
  attr(CR2_iv, "inverse_var") <- FALSE
  attr(CR2_iv, "target") <- attr(CR2_not, "target")
  expect_equal(CR2_not, CR2_iv)
  
  expect_identical(vcovCR(m3_hier, type = "CR2", target = iv), CR2_not)
  expect_identical(vcovCR(m3_hier, type = "CR2", target = iv, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(m3_hier, type = "CR2", target = m3_hier$data.full$var.eff.size), CR2_not))
})

