context("Wald tests")
set.seed(20190513)

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

lm_urbanicity <- lm(math1 ~ schoolk * stark + gender + ethnicity + lunchk, 
                    data = STAR)
V_urbanicity <- vcovCR(lm_urbanicity, cluster = STAR$schoolidk, type = "CR2")

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
