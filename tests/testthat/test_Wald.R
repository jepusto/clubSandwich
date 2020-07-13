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
  expect_identical(constrain_equal(constraints_cha)(coefs_sep), constraints_mat)
  
  
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
  
})

test_that("constrain_zero expressions are equivalent", {
  
  constraints_lgl <- grepl("typeprof:", coef_names_int)  
  constraints_int <- which(constraints_lgl)
  constraints_num <- as.numeric(constraints_int)
  constraints_char <- coef_names_int[constraints_lgl]
  constraints_mat <- diag(1L, nrow = length(coef_names_int))[constraints_lgl,,drop=FALSE]

  expect_identical(constrain_zero("typeprof:", coefs_int, reg_ex = TRUE), constraints_mat)
  expect_identical(constrain_zero(constraints_lgl, coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_int, coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_num, coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_char, coefs_int), constraints_mat)
  
  expect_type(constrain_zero("typeprof:", reg_ex = TRUE), "closure")
  expect_identical(constrain_zero("typeprof:", reg_ex = TRUE)(coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_lgl)(coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_int)(coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_num)(coefs_int), constraints_mat)
  expect_identical(constrain_zero(constraints_char)(coefs_int), constraints_mat)
  
})

test_that("constraint expressions are equivalent across specifications", {
  

})

test_that("Wald test is equivalent to Satterthwaite for q = 1.",{
  p <- length(coefs)
  t_tests <- coef_test(duncan_fit, vcov = Duncan_CR2)
  F_tests <- Wald_test(duncan_fit, vcov = Duncan_CR2, 
                       constraints = as.list(1:9))
  expect_equal(t_tests$df, sapply(F_tests, function(x) x$df))
  expect_equal(t_tests$p_Satt, sapply(F_tests, function(x) x$p_val))
})
