context("population mean estimation")
set.seed(20190513)

m <- 14
icc <- 0.2
mu <- 5
size <- 2
nj <- 1 + rnbinom(m, size = size, mu = mu)
group <- factor(rep(LETTERS[1:m], nj))
N <- sum(nj)

Y <- rnorm(m, sd = sqrt(icc))[group] + rnorm(N, sd = sqrt(1 - icc))
y_bar <- tapply(Y, group, mean)
lm_fit <- lm(Y ~ 1)

test_that("CR0 and df agree with formulas", {
  CR0 <- coef_test(lm_fit, vcov = "CR0", cluster = group, test = "Satterthwaite")
  VCR0_f <- sum(nj^2 * (y_bar - mean(Y))^2) / sum(nj)^2
  df0_f <- (N^2 - sum(nj^2))^2 / (N^2 * sum(nj^2) - 2 * N * sum(nj^3) + sum(nj^2)^2)
  
  expect_equal(as.numeric(CR0$SE), sqrt(VCR0_f))
  expect_equal(CR0$df, df0_f)
})

test_that("CR1 and df agree with formulas", {
  CR1 <- coef_test(lm_fit, vcov = "CR1", cluster = group, test = "Satterthwaite")
  VCR1_f <- (m / (m - 1)) * sum(nj^2 * (y_bar - mean(Y))^2) / sum(nj)^2
  df1_f <- (N^2 - sum(nj^2))^2 / (N^2 * sum(nj^2) - 2 * N * sum(nj^3) + sum(nj^2)^2)
  
  expect_equal(as.numeric(CR1$SE), sqrt(VCR1_f))
  expect_equal(CR1$df, df1_f)
})

test_that("CR2 and df agree with formulas", {
  CR2 <- coef_test(lm_fit, vcov = "CR2", cluster = group, test = "Satterthwaite")
  VCR2_f <- sum(nj^2 * (y_bar - mean(Y))^2 / (1 - nj / N)) / sum(nj)^2
  df2_f <- N^2 / (N^2 * sum(nj^2 / (N - nj)^2) - 2 * N * sum(nj^3 / (N - nj)^2) + sum(nj^2 / (N - nj))^2)
  
  expect_equal(as.numeric(CR2$SE), sqrt(VCR2_f))
  expect_equal(CR2$df, df2_f)
})

test_that("CR3 agrees with formula", {
  CR3 <- coef_test(lm_fit, vcov = "CR3", cluster = group, test = "Satterthwaite")
  VCR3_f <- sum(nj^2 * (y_bar - mean(Y))^2 / (1 - nj / N)^2) / sum(nj)^2
  # df2_f <- N^2 / (N^2 * sum(nj^2 / (N - nj)^2) - 2 * N * sum(nj^3 / (N - nj)^2) + sum(nj^2 / (N - nj))^2)
  
  expect_equal(as.numeric(CR3$SE), sqrt(VCR3_f))
  # expect_equal(CR2$df, df2_f)
})

test_that("CR4 and df agree with formulas", {
  CR4 <- coef_test(lm_fit, vcov = "CR4", cluster = group, test = "Satterthwaite")
  VCR4_f <- sum(nj^2 * (y_bar - mean(Y))^2 / (1 - nj / N)) / sum(nj)^2
  df4_f <- N^2 / (N^2 * sum(nj^2 / (N - nj)^2) - 2 * N * sum(nj^3 / (N - nj)^2) + sum(nj^2 / (N - nj))^2)
  
  expect_equal(as.numeric(CR4$SE), sqrt(VCR4_f))
  expect_equal(CR4$df, df4_f)
})
