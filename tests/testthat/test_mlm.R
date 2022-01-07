context("mlm objects")
set.seed(20190513)

n <- nrow(iris)
lm_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length, data = iris)
lm_A_fit <- lm(Sepal.Length ~ Species + Petal.Length, data = iris)
lm_B_fit <- lm(Sepal.Width ~ Species + Petal.Length, data = iris)
WLS_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length, data = iris, weights = Petal.Width)
CR_types <- paste0("CR",0:4)

test_that("bread works", {
  
  expect_equal(bread.mlm(lm_fit), sandwich:::bread.mlm(lm_fit))
  
  y <- with(iris, as.vector(rbind(Sepal.Length, Sepal.Width)))
  cluster <- rep(rownames(iris), each = ncol(residuals(lm_fit)))
  expect_true(check_bread(lm_fit, cluster = cluster, y = y))
  
  expect_true(check_bread(WLS_fit, cluster = cluster, y = y))
  
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(lm_fit, vcov = "CR2"))
  expect_true(check_CR(WLS_fit, vcov = "CR2"))
  expect_true(check_CR(lm_fit, vcov = "CR4"))
  expect_true(check_CR(WLS_fit, vcov = "CR4"))
})

test_that("vcovCR is mostly equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)
  CR_mats <- sapply(c("CR0","CR2","CR3","CR1","CR1p","CR1S"), 
                    function(t) as.matrix(vcovCR(lm_fit, type = t)),
                    simplify = FALSE, USE.NAMES = TRUE)
  HC_mats <- sapply(c("HC0","HC2","HC3","HC1"), 
                    function(t) vcovHC(lm_fit, type = t),
                    simplify = FALSE, USE.NAMES = TRUE)
  
  expect_equal(CR_mats$CR0, HC_mats$HC0)
  expect_equal(CR_mats$CR2, HC_mats$HC2)
  expect_equal(CR_mats$CR3, HC_mats$HC3)
  
  J <- nobs(lm_fit)
  p <- ncol(model.matrix(lm_fit))
  N <- nrow(model_matrix(lm_fit))
  expect_equal(CR_mats$CR1 * (J - 1), HC_mats$HC1 * (J - p))
  expect_equal(CR_mats$CR1p * (J - 2 * p), HC_mats$HC1 * (J - p))
  expect_equal(CR_mats$CR1S * (J - 1) * (N - 2 * p) / (N - 1), HC_mats$HC1 * (J - p))
  
  HC_A_mats <- sapply(c("HC0","HC2","HC3"), 
                      function(t) vcovHC(lm_A_fit, type = t),
                      simplify = FALSE, USE.NAMES = TRUE)
  HC_B_mats <- sapply(c("HC0","HC2","HC3"), 
                      function(t) vcovHC(lm_B_fit, type = t),
                      simplify = FALSE, USE.NAMES = TRUE)
  expect_equal(CR_mats$CR0[1:p,1:p], HC_A_mats$HC0, check.attributes = FALSE)
  expect_equal(CR_mats$CR2[1:p,1:p], HC_A_mats$HC2, check.attributes = FALSE)
  expect_equal(CR_mats$CR3[1:p,1:p], HC_A_mats$HC3, check.attributes = FALSE)
  expect_equal(CR_mats$CR0[p + 1:p,p + 1:p], HC_B_mats$HC0, check.attributes = FALSE)
  expect_equal(CR_mats$CR2[p + 1:p,p + 1:p], HC_B_mats$HC2, check.attributes = FALSE)
  expect_equal(CR_mats$CR3[p + 1:p,p + 1:p], HC_B_mats$HC3, check.attributes = FALSE)
})

test_that("mlm is equivalent to lm with long data.", {
  iris_long <- reshape(iris, c("Sepal.Length","Sepal.Width"), 
                       direction = "long", times = "outcome")
  iris_long$outcome <- paste0("Sepal.", iris_long$time)
  
  lm_long <- lm(Sepal ~ 0 + outcome +  outcome:Species + outcome:Petal.Length, data = iris_long)
  i <- order(rep(1:2, 4))
  expect_equal(coef_CS(lm_fit), coef(lm_long)[i], check.attributes = FALSE)
  
  CR_fit <- lapply(CR_types, function(x) as.matrix(vcovCR(lm_fit, type = x)))
  CR_long <- lapply(CR_types, function(x) vcovCR(lm_long, type = x, cluster = iris_long$id)[i,i])
  expect_equivalent(CR_fit, CR_long)
  
  test_fit <- lapply(CR_types, function(x) coef_test(lm_fit, vcov = x, test = "All", p_values = FALSE))
  test_long <- lapply(CR_types, function(x) coef_test(lm_long, vcov = x, cluster = iris_long$id, test = "All", p_values = FALSE)[i,])
  compare_ttests(test_fit, test_long)
  
  CR_fit <- lapply(CR_types, function(x) as.matrix(vcovCR(lm_fit, type = x, cluster = iris$Petal.Length)))
  CR_long <- lapply(CR_types, function(x) vcovCR(lm_long, type = x, cluster = iris_long$Petal.Length)[i,i])
  expect_equivalent(CR_fit, CR_long)
  
  test_fit <- lapply(CR_types, function(x) coef_test(lm_fit, vcov = x, test = "All", p_values = FALSE))
  test_long <- lapply(CR_types, function(x) coef_test(lm_long, vcov = x, cluster = iris_long$id, test = "All", p_values = FALSE)[i,])
  compare_ttests(test_fit, test_long)
  
})

test_that("Order doesn't matter.",{
  
  skip_on_cran()
  check_sort_order(WLS_fit, iris)
  
})

test_that("clubSandwich works with dropped covariates", {
  dat_miss <- iris
  dat_miss$Petal.Length[sample.int(n, size = round(n / 10))] <- NA
  lm_dropped <- update(lm_fit, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(Petal.Length))
  lm_complete <- update(lm_fit, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


test_that("clubSandwich works with dropped outcomes", {
  dat_miss <- iris
  n <- nrow(dat_miss)
  dat_miss$Sepal.Length[sample.int(n, size = round(n / 10))] <- NA
  dat_miss$Sepal.Width[sample.int(n, size = round(n / 10))] <- NA
  lm_dropped <- update(lm_fit, data = dat_miss)
  dat_complete <- subset(dat_miss, !is.na(Sepal.Length) & !is.na(Sepal.Width))
  lm_complete <- update(lm_fit, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, test = "All", p_values = FALSE))
  compare_ttests(test_drop, test_complete)
})

test_that("clubSandwich works with dropped outcomes, covariates, and weights", {
  dat_miss <- iris
  n <- nrow(dat_miss)
  dat_miss$Sepal.Length[sample.int(n, size = round(n / 5))] <- NA
  dat_miss$Sepal.Width[sample.int(n, size = round(n / 5))] <- NA
  dat_miss$Petal.Length[sample.int(n, size = round(n / 5))] <- NA
  dat_miss$Petal.Width[sample.int(n, size = round(n / 5))] <- NA
  WLS_dropped <- update(WLS_fit, data = dat_miss)
  dat_complete <- subset(dat_miss, 
                         !is.na(Petal.Length) & !is.na(Petal.Width) & 
                           !is.na(Sepal.Length) & !is.na(Sepal.Width))
  WLS_complete <- update(WLS_fit, data = dat_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(WLS_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(WLS_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(WLS_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(WLS_complete, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})

test_that("weight scale doesn't matter", {
  
  lm_fit_w <- update(lm_fit, weights = rep(10, nrow(iris)))
  
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit, type = x))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit_w, type = x))
  
  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix))  
  
  target <- rep(1 + rpois(nrow(iris), lambda = 8), each = ncol(residuals(lm_fit)))
  unweighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit, type = x, target = target))
  weighted_fit <- lapply(CR_types, function(x) vcovCR(lm_fit_w, type = x, target = target * 15))
  
  expect_equal(lapply(unweighted_fit, as.matrix), 
               lapply(weighted_fit, as.matrix))  
  
})

test_that("clubSandwich works with weights of zero.", {
  
  data("LifeCycleSavings")
  n_life <- nrow(LifeCycleSavings)
  
  LifeCycleSavings$cl <- substr(rownames(LifeCycleSavings), 1, 1)
  table(LifeCycleSavings$cl)
  LifeCycleSavings$wt <- rpois(nrow(LifeCycleSavings), lambda = 0.8)
  table(LifeCycleSavings$wt)
  
  lm_full <- lm(cbind(dpi, ddpi) ~ pop15 + pop75 + sr, data = LifeCycleSavings, weights = wt)
  LCS_sub <- subset(LifeCycleSavings, wt > 0)
  lm_sub <- lm(cbind(dpi, ddpi) ~ pop15 + pop75 + sr, data = LCS_sub, weights = wt)
  
  CR_full <- lapply(CR_types, function(x) vcovCR(lm_full, type = x))
  CR_sub <- lapply(CR_types, function(x) vcovCR(lm_sub, type = x))
  expect_equal(CR_full, CR_sub, check.attributes = FALSE)
  
  test_full <- lapply(CR_types, function(x) coef_test(lm_full, vcov = x, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  test_sub <- lapply(CR_types, function(x) coef_test(lm_sub, vcov = x, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  expect_equal(test_full, test_sub, check.attributes = FALSE)
  
  CR_full <- lapply(CR_types, function(x) vcovCR(lm_full, cluster = LifeCycleSavings$cl, type = x))
  CR_sub <- lapply(CR_types, function(x) vcovCR(lm_sub, cluster = LCS_sub$cl, type = x))
  expect_equal(CR_full, CR_sub, check.attributes = FALSE)
  
  test_full <- lapply(CR_types, function(x) coef_test(lm_full, vcov = x, cluster = LifeCycleSavings$cl, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  test_sub <- lapply(CR_types, function(x) coef_test(lm_sub, vcov = x, cluster = LCS_sub$cl, test = c("z","naive-t","Satterthwaite"), p_values = TRUE))
  expect_equal(test_full, test_sub, check.attributes = FALSE)
  
  dat_miss <- LifeCycleSavings
  miss_indicator <- sample.int(n_life, size = round(n_life / 5))
  dat_miss$pop15[miss_indicator] <- NA
  dat_miss$cl[miss_indicator] <- NA
  with(dat_miss, table(wt, is.na(pop15)))
  
  lm_dropped <- lm(cbind(dpi, ddpi) ~ pop15 + pop75 + sr, data = dat_miss, weights = wt)
  dat_complete <- subset(dat_miss, !is.na(pop15))
  lm_complete <- lm(cbind(dpi, ddpi) ~ pop15 + pop75 + sr, data = dat_complete, weights = wt)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(lm_dropped, cluster = dat_miss$cl, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(lm_complete, cluster = dat_complete$cl, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(lm_dropped, vcov = x, cluster = dat_miss$cl, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(lm_complete, vcov = x, cluster = dat_complete$cl, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
  
})
