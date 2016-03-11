context("robu objects")

library(robumeta)
data(corrdat)

test_that("CR0 z-tests agree with robumeta for correlated effects", {
  corr_large <- robu(effectsize ~ males + college + binge, data = corrdat, 
                         modelweights = "CORR", studynum = studyid,
                         var.eff.size = var, small = FALSE)
  p <- length(coef_CR(corr_large))
  N <- corr_large$N
  robu_CR0 <- vcovCR(corr_large, type = "CR0")
  ztests <- coef_test(corr_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(corr_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(corr_large$reg_table$SE, ztests$SE)
  expect_equal(with(corr_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
                    ztests$p_z)
})

test_that("CR2 t-tests agree with robumeta for correlated effects", {
  corr_small <- robu(effectsize ~ males + college + binge, data = corrdat, 
                         modelweights = "CORR", studynum = studyid,
                         var.eff.size = var)
  robu_CR2 <- vcovCR(corr_small, type = "CR2")
  CR2_ttests <- coef_test(corr_small, vcov = robu_CR2, test = "Satterthwaite")
    
  expect_equivalent(corr_small$VR.r, as.matrix(robu_CR2))
  expect_equal(corr_small$dfs, CR2_ttests$df)
  expect_equal(corr_small$reg_table$prob, CR2_ttests$p_Satt)
})

data(hierdat)

test_that("CR0 z-tests agree with robumeta for hierarchical effects", {
  hier_large <- robu(effectsize ~ binge + followup + sreport + age,
                         data = hierdat, studynum = studyid,
                         var.eff.size = var, modelweights = "HIER", small = FALSE)
  p <- length(coef_CR(hier_large))
  N <- hier_large$N
  robu_CR0 <- vcovCR(hier_large, type = "CR0")
  ztests <- coef_test(hier_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(hier_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(hier_large$reg_table$SE, ztests$SE)
  expect_equal(with(hier_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
               ztests$p_z)
})

test_that("CR2 t-tests agree with robumeta for hierarchical effects", {
  hier_small <- robu(effectsize ~ binge + followup + sreport + age,
                         data = hierdat, studynum = studyid,
                         var.eff.size = var, modelweights = "HIER")
  robu_CR2 <- vcovCR(hier_small, type = "CR2")
  CR2_ttests <- coef_test(hier_small, vcov = robu_CR2, test = "Satterthwaite")
  
  expect_equivalent(hier_small$VR.r, as.matrix(robu_CR2))
  expect_equal(hier_small$dfs, CR2_ttests$df)
  expect_equal(hier_small$reg_table$prob, CR2_ttests$p_Satt)
})

hierdat$user_wt <- 1 + rpois(nrow(hierdat), lambda = 3)

test_that("CR0 z-tests agree with robumeta for user weighting", {
  user_large <- robu(effectsize ~ binge + followup + sreport + age,
                     data = hierdat, studynum = studyid,
                     var.eff.size = var, userweights = user_wt, small = FALSE)
  p <- length(coef_CR(user_large))
  N <- user_large$N
  robu_CR0 <- vcovCR(user_large, type = "CR0")
  ztests <- coef_test(user_large, vcov = robu_CR0 * N / (N - p), test = "z")
  
  expect_equivalent(user_large$VR.r, as.matrix(robu_CR0))
  expect_equivalent(user_large$reg_table$SE, ztests$SE)
  expect_equal(with(user_large$reg_table, 2 * pnorm(abs(b.r / SE),lower.tail=FALSE)), 
               ztests$p_z)
})

test_that("CR2 t-tests agree with robumeta for user weighting", {
  user_small <- robu(effectsize ~ binge + followup + sreport + age,
                     data = hierdat, studynum = studyid,
                     var.eff.size = var, userweights = user_wt)
  user_lm <- lm(effectsize ~ binge + followup + sreport + age, data = hierdat,
                weights = user_wt)
  expect_equivalent(coef_CR(user_lm), coef(user_lm))
  
  robu_CR2 <- vcovCR(user_small, type = "CR2")
  expect_equivalent(user_small$VR.r, as.matrix(robu_CR2))
  lm_CR2 <- vcovCR(user_lm, cluster = hierdat$studyid, type = "CR2", target = user_small$data.full$avg.var.eff.size)
  expect_equivalent(robu_CR2, lm_CR2)
  
  CR2_ttests <- coef_test(user_small, vcov = robu_CR2, test = "Satterthwaite")
#   expect_equal(user_small$dfs, CR2_ttests$df)
#   expect_equal(user_small$reg_table$prob, CR2_ttests$p_Satt)
  lm_CR2_ttests <- coef_test(user_lm, vcov = "CR2", 
                             cluster = hierdat$studyid, 
                             target = user_small$data.full$avg.var.eff.size,
                             test = "Satterthwaite")
  expect_equivalent(CR2_ttests, lm_CR2_ttests)
})



data(dropoutPrevention)

test_that("dropoutPrevention tests replicate Tipton & Pustejovsky (2015) - full sample", {
  m3_hier <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                  + outcome + evaluator_independence
                  + male_pct + white_pct + average_age
                  + implementation_quality + program_site + duration + service_hrs, 
                  data = dropoutPrevention, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")
  m3_hier_CR2 <- vcovCR(m3_hier, cluster = dropoutPrevention$studyID, type = "CR2")
  CR2_ttests <- coef_test(m3_hier, vcov = m3_hier_CR2, test = "Satterthwaite")
  
  expect_equivalent(m3_hier$VR.r, as.matrix(m3_hier_CR2))
  expect_equal(m3_hier$dfs, CR2_ttests$df)
  expect_equal(m3_hier$reg_table$prob, CR2_ttests$p_Satt)

  contrast_list <- list("Study design" = 2:3, 
                        "Outcome measure" = 7:9,
                        "Evaluator independence" = 10:12,
                        "Implmentation quality" = 16:17,
                        "Program format" = 18:20)
  
  dropout_tests <- Wald_test(m3_hier, constraints = contrast_list, 
                             vcov = m3_hier_CR2, test = c("Naive-F","HTZ"))
  
  Fstat_club <- sapply(dropout_tests, function(x) x$F)
  Fstat_paper <- matrix(c(0.23, 0.22, 0.91, 0.84, 3.11, 2.78, 14.15, 13.78, 3.85, 3.65), nrow = 2)
  expect_equivalent(Fstat_paper, round(Fstat_club, 2))
  
  df_club <- sapply(dropout_tests, function(x) x$df[2])
  df_paper <- c(42.9, 21.5, 16.8, 36.9, 37.5)
  expect_equivalent(df_paper, round(df_club, 1))
})

test_that("dropoutPrevention tests replicate Tipton & Pustejovsky (2015) - reduced sample", {
  dp_subset <- subset(dropoutPrevention, big_study==TRUE)
  m3_hier <- robu(LOR1 ~ study_design + attrition + group_equivalence + adjusted
                  + outcome + evaluator_independence
                  + male_pct + white_pct + average_age
                  + implementation_quality + program_site + duration + service_hrs, 
                  data = dp_subset, studynum = studyID, var.eff.size = varLOR, modelweights = "HIER")
  m3_hier_CR2 <- vcovCR(m3_hier, cluster = dp_subset$studyID, type = "CR2")
  CR2_ttests <- coef_test(m3_hier, vcov = m3_hier_CR2, test = "Satterthwaite")
  
  expect_equivalent(m3_hier$VR.r, as.matrix(m3_hier_CR2))
  expect_equal(m3_hier$dfs, CR2_ttests$df)
  expect_equal(m3_hier$reg_table$prob, CR2_ttests$p_Satt)
  
  contrast_list <- list("Study design" = 2:3, 
                        "Outcome measure" = 7:9,
                        "Evaluator independence" = 10:11,
                        "Implmentation quality" = 15:16,
                        "Program format" = 17:19)
  
  dropout_tests <- Wald_test(m3_hier, constraints = contrast_list, 
                             vcov = "CR2", test = c("Naive-F","HTZ"))
  
  Fstat_club <- sapply(dropout_tests, function(x) x$F)
  Fstat_paper <- matrix(c(3.19, 2.93, 1.05, 0.84, 0.32, 0.26, 4.02, 3.69, 1.19, 0.98), nrow = 2)
  expect_equivalent(Fstat_paper, round(Fstat_club, 2))
  
  df_club <- sapply(dropout_tests, function(x) x$df[2])
  df_paper <- c(11.0, 7.7, 4.6, 11.0, 9.1)
  expect_equivalent(df_paper, round(df_club, 1))
})

# test user weighting
# test target matrix specification
# test inverse_var detection