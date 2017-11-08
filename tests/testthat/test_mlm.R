context("mlm objects")

lm_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length, data = iris)
WLS_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length, data = iris, weights = Petal.Width)
CR_types <- paste0("CR",0:4)

obj <- lm_fit

test_that("bread works", {
  y <- with(iris, as.vector(rbind(Sepal.Length, Sepal.Width)))
  cluster <- rep(rownames(iris), each = ncol(residuals(lm_fit)))
  expect_true(check_bread(lm_fit, cluster = cluster, y = y))
  expect_true(check_bread(WLS_fit, cluster = cluster, y = y))
  
})
