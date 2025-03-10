context("Wald tests")
set.seed(20190513)


skip_if_not_installed("carData")

data(Duncan, package = "carData")
Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
Duncan_int <- lm(prestige ~ type * (income + education), data=Duncan)
coefs_int <- coef(Duncan_int)
coef_names_int <- names(coefs_int)
Duncan_int_CR2 <- vcovCR(Duncan_int, type = "CR2", cluster = Duncan$cluster)

Duncan_sep <- lm(prestige ~ 0 + type + type:(income + education), data=Duncan)
coefs_sep <- coef(Duncan_sep)
coef_names_sep <- names(coefs_sep)
Duncan_sep_CR2 <- vcovCR(Duncan_sep, type = "CR2", cluster = Duncan$cluster)

test_that("constrain_equal expressions are equivalent", {
  
  constraints_lgl <- grepl(":education", coef_names_sep)  
  constraints_int <- which(constraints_lgl)
  constraints_num <- as.numeric(constraints_int)
  constraints_char <- coef_names_sep[constraints_lgl]
  constraints_mat <- cbind(matrix(0L, 2, 6), matrix(c(-1L, -1L, 1L, 0L, 0L, 1L), 2, 3))
  
  expect_identical(constrain_equal(":education", coefs_sep, reg_ex = TRUE), constraints_mat)
  expect_identical(constrain_equal(constraints_lgl, coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_int, coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_num, coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_char, coefs_sep), constraints_mat)

  expect_type(constrain_equal(":education", reg_ex = TRUE), "closure")
  expect_identical(constrain_equal(":education", reg_ex = TRUE)(coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_lgl)(coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_int)(coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_num)(coefs_sep), constraints_mat)
  expect_identical(constrain_equal(constraints_char)(coefs_sep), constraints_mat)
  
  constraint_list <- constrain_equal(list(type = 1:3, income = 4:6, edu = 7:9),
                                    coefs = coefs_sep) 
  constraint_func <- constrain_equal(list(type = 1:3, income = 4:6, edu = 7:9))
  expect_identical(constraint_list, constraint_func(coefs_sep))

  Wald_A <- Wald_test(Duncan_sep, constraints = constraint_list,
                      vcov = Duncan_sep_CR2, type = "All")
  Wald_B <- Wald_test(Duncan_sep, constraints = constraint_func,
                      vcov = Duncan_sep_CR2, type = "All")
  expect_identical(Wald_A, Wald_B)    
})

test_that("constrain_pairwise expressions are equivalent", {
  
  constraints_lgl <- grepl(":education", coef_names_sep)  
  constraints_int <- which(constraints_lgl)
  constraints_num <- as.numeric(constraints_int)
  constraints_char <- coef_names_sep[constraints_lgl]
  constraints_mat <- constrain_pairwise(":education", coefs_sep, reg_ex = TRUE)
  
  expect_identical(length(constraints_mat), sum(constraints_lgl))
  expect_identical(constrain_pairwise(constraints_lgl, coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_int, coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_num, coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_char, coefs_sep), constraints_mat)
  
  expect_type(constrain_pairwise(":education", reg_ex = TRUE), "closure")
  expect_identical(constrain_pairwise(constraints_lgl)(coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_int)(coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_num)(coefs_sep), constraints_mat)
  expect_identical(constrain_pairwise(constraints_char)(coefs_sep), constraints_mat)
  
  constraint_list <- constrain_pairwise(list(type = 1:3, income = 4:6, edu = 7:9),
                                     coefs = coefs_sep) 
  constraint_func <- constrain_pairwise(list(type = 1:3, income = 4:6, edu = 7:9))
  expect_identical(constraint_list, constraint_func(coefs_sep))
  
  Wald_A <- Wald_test(Duncan_sep, constraints = constraint_list,
                      vcov = Duncan_sep_CR2, type = "All")
  Wald_B <- Wald_test(Duncan_sep, constraints = constraint_func,
                      vcov = Duncan_sep_CR2, type = "All")
  expect_identical(Wald_A, Wald_B)
})

test_that("constrain_zero expressions are equivalent", {
  
  constraints_lgl <- grepl("typeprof:", coef_names_int)  
  constraints_int <- which(constraints_lgl)
  constraints_num <- as.numeric(constraints_int)
  constraints_char <- coef_names_int[constraints_lgl]
  constraints_mat <- diag(1L, nrow = length(coef_names_int))[constraints_lgl,,drop=FALSE]

  expect_equal(constrain_zero("typeprof:", coefs_int, reg_ex = TRUE), constraints_mat)
  expect_equal(constrain_zero(constraints_lgl, coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_int, coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_num, coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_char, coefs_int), constraints_mat)
  
  expect_type(constrain_zero("typeprof:", reg_ex = TRUE), "closure")
  expect_equal(constrain_zero("typeprof:", reg_ex = TRUE)(coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_lgl)(coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_int)(coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_num)(coefs_int), constraints_mat)
  expect_equal(constrain_zero(constraints_char)(coefs_int), constraints_mat)
  
  constraint_list <- constrain_zero(list(type = 2:3, income = 6:7, edu = 8:9),
                                        coefs = coefs_int) 
  constraint_func <- constrain_zero(list(type = 2:3, income = 6:7, edu = 8:9))
  expect_equal(constraint_list, constraint_func(coefs_int))
  
  Wald_A <- Wald_test(Duncan_int, constraints = constraint_list,
                      vcov = Duncan_int_CR2, type = "All")
  Wald_B <- Wald_test(Duncan_int, constraints = constraint_func,
                      vcov = Duncan_int_CR2, type = "All")
  expect_equal(Wald_A, Wald_B)
})

test_that("constraint expressions are equivalent across specifications", {

  skip_on_cran()
  skip_if(R.version$major < "4", "Skip for R versions below 4.")
  
  constraints_eq <- constrain_equal(
    list(type = 1:3, income = 4:6, edu = 7:9),
    coefs = coefs_sep
  )  
  # constraints_eq$all <- do.call(rbind, constraints_eq)
  
  constraints_null <- constrain_zero(
    list(type = 2:3, income = 6:7, edu = 8:9),
    coefs = coefs_int
  )
  # constraints_null$all <- do.call(rbind, constraints_null)
  
  Wald_eq <- Wald_test(Duncan_sep, 
                       constraints_eq, 
                       vcov = Duncan_sep_CR2, 
                       test = c("Naive-F","HTZ","EDF"), tidy = TRUE) 

  Wald_zero <- Wald_test(Duncan_int, 
                         constraints_null,
                         vcov = Duncan_int_CR2, 
                         test = c("Naive-F","HTZ","EDF"), tidy = TRUE)
  
  compare_Waldtests(Wald_eq, Wald_zero)

  pairwise_sep <- constrain_pairwise(
    list(type = 1:3, income = 4:6, edu = 7:9),
    coefs = coefs_sep
  )
  
  pairwise_int <- constrain_pairwise(
    list(type = 2:3, income = 6:7, edu = 8:9),
    coefs = coefs_int,
    with_zero = TRUE
  )
  
  pairwise_sep <- Wald_test(Duncan_sep, 
                       pairwise_sep, 
                       vcov = Duncan_sep_CR2, 
                       tidy = TRUE) 
  
  pairwise_int <- Wald_test(Duncan_int, 
                         pairwise_int,
                         vcov = Duncan_int_CR2,
                         tidy = TRUE)
  
  compare_Waldtests(pairwise_sep, pairwise_int)
  
})

test_that("Wald test is equivalent to Satterthwaite for q = 1.", {
  
  skip_on_cran()
  
  t_tests_sep <- coef_test(Duncan_sep, vcov = Duncan_sep_CR2)
  
  constraints_sep <- as.list(1:9)
  names(constraints_sep) <- coef_names_sep
  
  F_tests_sep <- Wald_test(Duncan_sep, vcov = Duncan_sep_CR2, 
                           constraints = constrain_zero(constraints_sep),
                           tidy = TRUE)
      
  expect_equal(t_tests_sep$tstat^2, F_tests_sep$Fstat, tol = 10^-5)
  expect_equal(rep(1, 9), F_tests_sep$df_num, tol = 10^-5)
  expect_equal(t_tests_sep$df, F_tests_sep$df_denom, tol = 10^-5)
  expect_equal(t_tests_sep$p_Satt, F_tests_sep$p_val, tol = 10^-5)
  
  t_tests_int <- coef_test(Duncan_int, vcov = Duncan_int_CR2)
  
  constraints_int <- as.list(1:9)
  names(constraints_int) <- coef_names_int
  
  F_tests_int <- Wald_test(Duncan_int, vcov = Duncan_int_CR2, 
                           constraints = constrain_zero(constraints_int),
                           tidy = TRUE)
  
  expect_equal(t_tests_int$tstat^2, F_tests_int$Fstat, tol = 10^-5)
  expect_equal(rep(1, 9), F_tests_int$df_num, tol = 10^-5)
  expect_equal(t_tests_int$df, F_tests_int$df_denom, tol = 10^-5)
  expect_equal(t_tests_int$p_Satt, F_tests_int$p_val, tol = 10^-5)
  
})


skip_if_not_installed("AER")
data(STAR, package = "AER")

# clean up a few variables
levels(STAR$stark)[3] <- "aide"
levels(STAR$schoolk)[1] <- "urban"
STAR <- subset(STAR, 
               !is.na(schoolidk),
               select = c(schoolidk, schoolk, stark, gender, ethnicity, math1, lunchk))

lm_urbanicity <- lm(math1 ~ schoolk * stark + gender + ethnicity + lunchk, 
                    data = STAR)
V_urbanicity <- vcovCR(lm_urbanicity, cluster = STAR$schoolidk, type = "CR2")

test_that("Wald_test works with lists.", {
  test_A <- Wald_test(lm_urbanicity, 
                      constraints = constrain_zero("schoolk.+:stark", reg_ex = TRUE),
                      vcov = V_urbanicity)
  
  test_B <- Wald_test(lm_urbanicity, 
                      constraints = constrain_zero("schoolk.+:starksmall", reg_ex = TRUE),
                      vcov = V_urbanicity)
  
  C_list <- list(
    `Any interaction` = constrain_zero("schoolk.+:stark", 
                                       coef(lm_urbanicity), reg_ex = TRUE),
    `Small vs regular` = constrain_zero("schoolk.+:starksmall", 
                                        coef(lm_urbanicity), reg_ex = TRUE)
  )
  
  D_list <- constrain_zero(constraints = list(
    `Any interaction` = "schoolk.+:stark",
    `Small vs regular` = "schoolk.+:starksmall"
    ), reg_ex = TRUE)
  
  test_C <- Wald_test(lm_urbanicity, 
                          constraints = C_list,
                          vcov = V_urbanicity)
  
  test_D <- Wald_test(lm_urbanicity, 
                      constraints = D_list,
                      vcov = V_urbanicity)
  
  test_E <- Wald_test(
                    lm_urbanicity, 
                    constraints = list(
                      `Any interaction` = constrain_zero("schoolk.+:stark", reg_ex = TRUE),
                      `Small vs regular` = constrain_zero("schoolk.+:starksmall", reg_ex = TRUE)
                    ),
                    vcov = V_urbanicity
                  )
  
  expect_identical(test_A, test_C$`Any interaction`)
  expect_identical(test_A, test_D$`Any interaction`)
  expect_identical(test_A, test_E$`Any interaction`)
  expect_identical(test_B, test_C$`Small vs regular`)
  expect_identical(test_B, test_D$`Small vs regular`)
  expect_identical(test_B, test_E$`Small vs regular`)
  
})

test_that("Wald_test has informative error messages.", {
  
  expect_error(
    Wald_test(lm_urbanicity, 
              constraints = constrain_zero("schoolk.+:stark", reg_ex = TRUE),
              vcov = V_urbanicity, 
              test = "none"
              )
  )
  
  A <- Wald_test(lm_urbanicity, 
            constraints = constrain_zero("schoolk.+:stark", reg_ex = TRUE),
            vcov = V_urbanicity, 
            test = c("none","HTA")
  )

  B <- Wald_test(lm_urbanicity, 
                 constraints = constrain_zero("schoolk.+:stark", reg_ex = TRUE),
                 vcov = V_urbanicity, 
                 test = "All"
  )
  
  expect_equal(A, subset(B, test == "HTA"), check.attributes = FALSE)
})


test_that("Wald_test works for intercept-only models.", {
  
  lm_int <- lm(math1 ~ 1, data = STAR)
  vcov_int <- vcovCR(lm_int, cluster = STAR$schoolidk, type = "CR2")
  F_test <- Wald_test(lm_int, constraints = constrain_zero(1), 
                      vcov = vcov_int, test = c("HTZ","HTA","HTB"))
  t_test <- coef_test(lm_int, vcov = vcov_int)
  
  expect_equal(F_test$Fstat, rep(t_test$tstat^2, 3L))
  expect_equal(F_test$df_denom, rep(t_test$df, 3L))
  expect_equal(F_test$p_val, rep(t_test$p_Satt, 3L))
  
  lm_sep <- lm(math1 ~ 0 + schoolk, data = STAR)
  vcov_sep <- vcovCR(lm_sep, cluster = STAR$schoolidk, type = "CR2")
  F_test <- Wald_test(lm_sep, 
                      constraints = constrain_pairwise(1:3, with_zero = TRUE), 
                      vcov = vcov_sep, test = "HTZ", tidy = TRUE)
  
  t_test <- coef_test(lm_sep, vcov = vcov_sep)
  
  expect_equal(F_test$Fstat[1:3], t_test$tstat^2)
  expect_equal(F_test$df_denom[1:3], t_test$df)
  expect_equal(F_test$p_val[1:3], t_test$p_Satt)
  
})

test_that("Wald_test fails gracefully when between-cluster variance of coefficients isn't identified.", {

  skip_if_not_installed("metafor")
  
  suppressPackageStartupMessages(library(metafor))
  
  dat <- dat.bcg
  dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat, subset=-5)
  res <- rma(yi, vi, data=dat, mods = ~ 0 + alloc)
  
  Vmat <- vcovCR(res, cluster=dat$trial, type="CR2")
  expect_equal(Vmat[1,1], 0)
  
  t_tests <- coef_test(res, cluster=dat$trial, vcov="CR2")
  expect_true(is.na(t_tests$df_Satt[1]))
  expect_true(is.na(t_tests$p_Satt[1]))
  
  CI <- conf_int(res, cluster=dat$trial, vcov="CR2")
  expect_true(is.na(CI$CI_L[1]))
  expect_true(is.na(CI$CI_U[1]))
  
  Wald1 <- Wald_test(
    res, 
    cluster=dat$trial, 
    vcov="CR2", 
    constraints=constrain_equal(1:3),
    test = "All"
  )
  
  expect_s3_class(Wald1, "Wald_test_clubSandwich")
  
  expect_error(
    Wald_test(
      res, 
      cluster=dat$trial, 
      vcov="CR2", 
      constraints=constrain_zero(1:3)
    ), 
    regexp = "not positive definite"
  )
  
  Wald2 <- Wald_test(
    res, 
    cluster=dat$trial, 
    vcov="CR2", 
    constraints = list(A = constrain_equal(1:3), B = constrain_zero(1:3)),
    test = "All"
  )
  
  expect_s3_class(Wald2$A, "Wald_test_clubSandwich") 
  expect_s3_class(Wald2$B, "Wald_test_clubSandwich") 
  expect_identical(Wald1, Wald2$A)  
  expect_true(all(is.na(Wald2$B$Fstat)))
  expect_true(all(is.na(Wald2$B$p_val)))
  
  Wald3 <- Wald_test(
    res, 
    cluster=dat$trial, 
    vcov="CR2", 
    constraints = list(A = constrain_equal(1:3), B = constrain_zero(1:3)),
    tidy = TRUE
  )
  
  expect_s3_class(Wald3, "Wald_test_clubSandwich") 
  expect_equivalent(Wald1[Wald1$test=="HTZ",], Wald3[1,-1])  
  expect_true(is.na(Wald3[2,"Fstat"]))
  expect_true(is.na(Wald3[2,"p_val"]))
  
  Wald4 <- Wald_test(
    res, 
    cluster=dat$trial, 
    vcov="CR2", 
    constraints=constrain_pairwise(1:3, with_zero = TRUE),
    tidy = TRUE
  )
  
  expect_s3_class(Wald4, "Wald_test_clubSandwich")
  expect_true(is.na(Wald4[1,"Fstat"]))
  expect_true(is.na(Wald4[1,"p_val"]))
  
})

test_that("Wald_test words with multiple comparisons adjustment", {
  
  Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
    
  Wald5 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTZ","chi-sq")
  )
  
  Wald6 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTZ","chi-sq"),
    adjustment_method = "none"
  )

  # check that explicitly stating default does not affect functionality
  lapply(Wald5, expect_s3_class, class = "Wald_test_clubSandwich")
  lapply(Wald6, expect_s3_class, class = "Wald_test_clubSandwich")
  expect_equal(Wald5, Wald6) 

  # change Wald6 to have hochberg adjustment
  Wald6 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTZ","chi-sq"),
    adjustment_method = "hochberg"
  )

  Wald5_p_values <- sapply(Wald5, function(x) x$p_val) # extract p-values of Wald5
  Wald6_p_values <- sapply(Wald6, function(x) x$p_val) # extract p-values of Wald6
  
  expect_false(all(Wald5_p_values == Wald6_p_values))
  
  # get adjusted_p_values for Wald5, formatted the same as extracted p-values from Wald6
  Wald5_adjusted_p <- apply(Wald5_p_values, 1, p.adjust, method = "hochberg", simplify = FALSE)
  Wald5_adjusted_p <- do.call(rbind, Wald5_adjusted_p)
  expect_equal(Wald5_adjusted_p, Wald6_p_values)
  
  # Now using tidy = TRUE
  Wald7 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTA","EDF","EDT"),
    tidy = TRUE
  )
  
  Wald8 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTA","EDF","EDT"),
    adjustment_method = "none",
    tidy = TRUE
  )
  
  expect_equal(Wald7,Wald8)
  
  Wald9 <- Wald_test(
    Duncan_fit,
    constraints = constrain_pairwise(":education", reg_ex = TRUE),
    vcov = "CR2",
    cluster = Duncan$cluster,
    test = c("HTA","EDF","EDT"),
    adjustment_method = "holm",
    tidy = TRUE
  )
  
  Wald8_p_adjusted <- tapply(Wald8$p_val, Wald8$test, p.adjust, method = "holm", simplify = FALSE)
  Wald8_adjusted <- Wald8
  Wald8_adjusted$p_val <- unsplit(Wald8_p_adjusted, Wald8$test)
  
  expect_equal(Wald8_adjusted, Wald9)
  
})
