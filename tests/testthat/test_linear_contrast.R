context("linear contrasts")
set.seed(20210110)

skip_if_not_installed("carData")


# Duncan example
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


# STAR example
skip_if_not_installed("AER")
data(STAR, package = "AER")
levels(STAR$stark)[3] <- "aide"
levels(STAR$schoolk)[1] <- "urban"
STAR <- subset(STAR, 
               !is.na(schoolidk),
               select = c(schoolidk, schoolk, stark, gender, ethnicity, math1, lunchk))

lm_urbanicity <- lm(math1 ~ schoolk * stark + gender + ethnicity + lunchk, 
                    data = STAR)


CRs <- paste0("CR", 0:4)

test_that("vcov arguments work", {
  VCR <- lapply(CRs, function(t) vcovCR(Duncan_sep, type = t, cluster = Duncan$cluster))
  CI_A <- lapply(VCR, function(v) linear_contrast(Duncan_sep, 
                                                  vcov = v, cluster = Duncan$cluster, 
                                                  contrast = constrain_pairwise(1:3), 
                                                  level = .98))
  CI_B <- lapply(CRs, function(t) linear_contrast(Duncan_sep, 
                                                  vcov = t, cluster = Duncan$cluster, 
                                                  contrast = constrain_pairwise(1:3), 
                                                  level = .98))
  expect_equal(CI_A, CI_B, check.attributes = FALSE)
  
})

test_that("constrain_() functions work.", {
  
  # Not worrying about CR4 here
  
  CI_A <- lapply(CRs[1:4], function(t) 
    linear_contrast(Duncan_sep, vcov = t, cluster = Duncan$cluster,
                    contrast = constrain_pairwise(1:3))
  )
  
  CI_B <- lapply(CRs[1:4], function(t) 
    linear_contrast(Duncan_int, vcov = t, cluster = Duncan$cluster,
                    contrast = constrain_pairwise(2:3, with_zero = TRUE))
  )
  
  CI_A <- lapply(CI_A, subset, select = -Coef)
  CI_B <- lapply(CI_B, subset, select = -Coef)
  
  expect_equal(CI_A, CI_B, check.attributes = FALSE)
})

test_that("linear_contrast() works with lists.", {
  
  
  CIs <- linear_contrast(lm_urbanicity, 
                         vcov = "CR2", cluster = STAR$schoolidk,
                         contrast = list(
                           A = constrain_zero(2:3),
                           B = constrain_pairwise(2:3, with_zero = TRUE),
                           C = constrain_equal("ethnicity", reg_ex = TRUE),
                           D = constrain_pairwise("ethnicity", reg_ex = TRUE)
                          ))

  CI_A <- as.data.frame(subset(CIs, grepl("A\\.", Coef), select = -Coef))
  CI_B <- as.data.frame(subset(CIs, grepl("B\\.[a-z]+$", Coef), select = -Coef))
  expect_equal(CI_A, CI_B, check.attributes = FALSE)

  CI_C <- as.data.frame(subset(CIs, grepl("C\\.", Coef), select = -Coef))
  CI_D <- as.data.frame(subset(CIs, grepl("D\\..+ethnicityafam$", Coef), select = -Coef))
  expect_equal(CI_C, CI_D, check.attributes = FALSE)
  
})

test_that("printing works", {
  CIs <- linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                         constrain_pairwise(":education", reg_ex = TRUE))
  expect_output(print(CIs))
  
  CIs <- linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                         constrain_pairwise(":education", reg_ex = TRUE),
                         p_values = TRUE)
  expect_output(x <- print(CIs))
  expect_true(all(c("p-value","Sig.") %in% names(x)))
})

test_that("level checks work", {
  expect_error(linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                               constrain_pairwise(":education", reg_ex = TRUE),
                               level = -0.01))
  expect_error(linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                               constrain_pairwise(":education", reg_ex = TRUE),
                               level = 95))
  expect_output(print(
    linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                    constrain_pairwise(":education", reg_ex = TRUE),
                    level = runif(1))
  ))
})

test_that("CI boundaries are ordered", {
  lev <- runif(1)
  CI_z <- linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                          constrain_pairwise(":education", reg_ex = TRUE),
                          test = "z", level = lev)
  CI_t <- linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                          constrain_pairwise(":education", reg_ex = TRUE),
                          test = "naive-t", level = lev)
  CI_Satt <- linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                             constrain_pairwise(":education", reg_ex = TRUE),
                              test = "Satterthwaite", level = lev)
  expect_true(all(CI_t$CI_L < CI_z$CI_L))
  expect_true(all(CI_t$CI_U > CI_z$CI_U))
  expect_true(all(CI_Satt$CI_L < CI_z$CI_L))
  expect_true(all(CI_Satt$CI_U > CI_z$CI_U))
})

test_that("linear_contrast() is consistent with Wald_test()", {
  
  skip_on_cran()
  
  lev <- runif(1)
  CIs <- lapply(CRs, function(v) 
    linear_contrast(lm_urbanicity, vcov = v, cluster = STAR$schoolidk,
                    contrasts = list(
                      A = constrain_zero("ethnicity", reg_ex = TRUE),
                      B = constrain_equal("ethnicity", reg_ex = TRUE)
                    ),
                    test = "Satterthwaite", level = lev, p_values = TRUE))
  Wald_tests <- lapply(CRs, function(v) 
    Wald_test(lm_urbanicity, vcov = v, cluster = STAR$schoolidk,
              constraints = constrain_pairwise("ethnicity", reg_ex = TRUE, with_zero = TRUE),
              test = "HTZ", tidy = TRUE))
  
  CI_pvals <- lapply(CIs, function(x) x$p_val)
  Wald_pvals <- lapply(Wald_tests, function(x) x$p_val[1:7])
  expect_equal(CI_pvals, Wald_pvals, tolerance = 1e-6)

})

test_that("linear_contrast() has informative error messages.", {
  expect_error(
    linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                    constrain_pairwise(":education", reg_ex = TRUE),
                    test = "all")
  )
  
  expect_error(
    linear_contrast(Duncan_sep, vcov = "CR0", cluster = Duncan$cluster,
                    constrain_pairwise(":education", reg_ex = TRUE),
                    test = "saddlepoint")
  )
  
})
