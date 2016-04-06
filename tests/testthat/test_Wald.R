context("Wald tests")

data(Duncan, package = "car")
Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
duncan_fit <- lm(prestige ~ type * (income + education), data=Duncan)
coefs <- names(coef(duncan_fit))
Duncan_CR2 <- vcovCR(duncan_fit, type = "CR2", cluster = Duncan$cluster)

test_that("constraint expressions are equivalent", {
  constraints_logical <- grepl("typeprof:", coefs)  
  constraints_int <- which(constraints_logical)
  constraints_num <- as.numeric(constraints_int)
  constraints_char <- coefs[constraints_logical]
  constraints_mat <- diag(1L, nrow = length(coefs))[constraints_logical,,drop=FALSE]
  Wald_logical <- Wald_test(duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
                            constraints = constraints_logical, test = "All")
  expect_output(Wald_logical, "")
  
  constraint_list <- list(integer = constraints_int,
                          numeric = constraints_num, 
                          char = constraints_char,
                          matrix = constraints_mat)
  Walds <- Wald_test(duncan_fit, vcov = "CR2", cluster = Duncan$cluster, 
                     constraints = constraint_list, test = "All")
  expect_identical(Wald_logical, Walds$integer)
  expect_identical(Wald_logical, Walds$numeric)
  expect_identical(Wald_logical, Walds$char)
  expect_identical(Wald_logical, Walds$matrix)
})

test_that("Wald test is equivalent to Satterthwaite for q = 1.",{
  p <- length(coefs)
  t_tests <- coef_test(duncan_fit, vcov = Duncan_CR2)
  F_tests <- Wald_test(duncan_fit, vcov = Duncan_CR2, 
                       constraints = as.list(1:9))
  expect_equal(t_tests$df, sapply(F_tests, function(x) x$df))
  expect_equal(t_tests$p_Satt, sapply(F_tests, function(x) x$p_val))
})
