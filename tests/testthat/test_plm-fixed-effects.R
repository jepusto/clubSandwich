context("plm objects - fixed effects")

library(plm)

data("Produc", package = "plm")
Produc$cluster <- sample(LETTERS[1:10], size = nrow(Produc), replace=TRUE)
Produc_scramble <- Produc[sample(nrow(Produc)),]

plm_individual <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                      data = Produc_scramble, index = c("state","year"), 
                      effect = "individual", model = "within")

test_that("individual effects agree with lm", {
  lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
  coef_names <- names(coef(plm_individual)) 
  coef_index <- names(coef(lm_individual)) %in% coef_names
  
  expect_equal(vcovCR(plm_individual, type="CR0")[coef_names,coef_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR0")[coef_index,coef_index])
  expect_equal(vcovCR(plm_individual, type="CR1")[coef_names,coef_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR1")[coef_index,coef_index])
  expect_equal(vcovCR(plm_individual, type="CR2")[coef_names,coef_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[coef_index,coef_index])
})

plm_time <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc_scramble, index = c("state","year"), 
                effect = "time", model = "within")

test_that("time effects agree with lm", {
  
  lm_time <- lm(log(gsp) ~ 0 + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
  coef_names <- names(coef(plm_time)) 
  coef_index <- names(coef(lm_time)) %in% coef_names
  
  expect_equal(vcovCR(plm_time, type="CR0")[coef_names,coef_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR0")[coef_index,coef_index])
  expect_equal(vcovCR(plm_time, type="CR1")[coef_names,coef_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR1")[coef_index,coef_index])
  expect_equal(vcovCR(plm_time, type="CR2")[coef_names,coef_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR2")[coef_index,coef_index])
})

plm_twoways <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                   data = Produc_scramble, index = c("state","year"), 
                   effect = "twoways", model = "within")

test_that("two-way effects agree with lm", {
  
  lm_twoways <- lm(log(gsp) ~ 0 + state + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
  coef_names <- names(coef(plm_twoways)) 
  coef_index <- names(coef(lm_twoways)) %in% coef_names
  
  # clustering on individual
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR0")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR0")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR1")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR1")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR2")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR2")[coef_index,coef_index])

  # clustering on time
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR0")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR0")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR1")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR1")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR2")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR2")[coef_index,coef_index])

  # clustering on a randomly generated factor
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR0")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR0")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR1")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR1")[coef_index,coef_index])
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR2")[coef_names,coef_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR2")[coef_index,coef_index])
  
})

plm_fit <- plm_individual
n <- nrow(Produc_scramble)

obj <- plm_fit
type <- "CR2"
inverse_var <- FALSE
target <- NULL
index <- attr(model.frame(obj),"index")
cluster <- switch(obj$args$effect,
                  individual = index[[1]],
                  time = index[[2]])

test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(plm_fit, type = "CR2")
  expect_identical(vcovCR(plm_fit, cluster = Produc_scramble$state, type = "CR2"), CR2_iv)
  expect_identical(vcovCR(plm_fit, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_identical(vcovCR(plm_fit, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_fit, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_identical(vcovCR(plm_fit, cluster = Produc_scramble$state, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_identical(vcovCR(plm_fit, type = "CR2", target = rep(1, n)), CR2_not)
  expect_identical(vcovCR(plm_fit, type = "CR2", target = rep(1, n), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(plm_fit, type = "CR2", target = 1 / Produc_scramble$emp), CR2_not))
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_fit, type = "CR4")
  expect_identical(vcovCR(plm_fit, cluster = Produc_scramble$state, type = "CR4"), CR4_iv)
  expect_identical(vcovCR(plm_fit, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_identical(vcovCR(plm_fit, type = "CR4", target = rep(1, n), inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_fit, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_identical(vcovCR(plm_fit, cluster = Produc_scramble$state, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_identical(vcovCR(plm_fit, type = "CR4", target = rep(1, n)), CR4_not)
  expect_identical(vcovCR(plm_fit, type = "CR4", target = rep(1, n), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(plm_fit, type = "CR4", target = 1 / Produc_scramble$emp), CR4_not))
  
})

test_that("CR2 and CR4 are target-unbiased", {
  
  expect_true(check_CR(plm_individual, vcov = "CR2"))
  expect_true(check_CR(plm_individual, vcov = "CR4"))
  
  expect_true(check_CR(plm_time, vcov = "CR2"))
  expect_true(check_CR(plm_time, vcov = "CR4"))
  
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "individual"))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "individual"))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "time"))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "time"))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = Produc_scramble$cluster))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = Produc_scramble$cluster))
})

test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich)
  CR0 <- vcovCR(plm_fit, cluster = 1:n, type = "CR0")
  expect_equal(vcovHC(plm_fit, type = "HC0"), as.matrix(CR0))
  CR1 <- vcovCR(plm_fit, cluster = 1:n, type = "CR1S")
  expect_equal(vcovHC(plm_fit, type = "HC1"), as.matrix(CR1) * (n - 1) / n)
  CR2 <- vcovCR(plm_fit, cluster = 1:n, type = "CR2")
  expect_equal(vcovHC(plm_fit, type = "HC2"), as.matrix(CR2))
  CR3 <- vcovCR(plm_fit, cluster = 1:n, type = "CR3")
  expect_equal(vcovHC(plm_fit, type = "HC3"), as.matrix(CR3))
})
