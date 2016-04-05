context("plm objects - fixed effects")

library(plm)

data("Produc", package = "plm")
Produc$cluster <- sample(LETTERS[1:10], size = nrow(Produc), replace=TRUE)
Produc_scramble <- Produc[sample(nrow(Produc)),]

plm_individual <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                      data = Produc_scramble, index = c("state","year"), 
                      effect = "individual", model = "within")
lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
individual_names <- names(coef(plm_individual)) 
individual_index <- names(coef(lm_individual)) %in% individual_names

test_that("individual effects agree with lm", {
  expect_equal(vcovCR(plm_individual, type="CR0")[individual_names,individual_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR0")[individual_index,individual_index])
  expect_equal(vcovCR(plm_individual, type="CR1")[individual_names,individual_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR1")[individual_index,individual_index])
  expect_equal(vcovCR(plm_individual, type="CR2")[individual_names,individual_names], 
               vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index])
})

plm_time <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                data = Produc_scramble, index = c("state","year"), 
                effect = "time", model = "within")
lm_time <- lm(log(gsp) ~ 0 + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
time_names <- names(coef(plm_time)) 
time_index <- names(coef(lm_time)) %in% time_names

test_that("time effects agree with lm", {
  expect_equal(vcovCR(plm_time, type="CR0")[time_names,time_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR0")[time_index,time_index])
  expect_equal(vcovCR(plm_time, type="CR1")[time_names,time_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR1")[time_index,time_index])
  expect_equal(vcovCR(plm_time, type="CR2")[time_names,time_names], 
               vcovCR(lm_time, cluster = Produc$year, type = "CR2")[time_index,time_index])
})

plm_twoways <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                   data = Produc_scramble, index = c("state","year"), 
                   effect = "twoways", model = "within")
lm_twoways <- lm(log(gsp) ~ 0 + state + factor(year) + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
twoway_names <- names(coef(plm_twoways)) 
twoway_index <- names(coef(lm_twoways)) %in% twoway_names


test_that("two-way effects agree with lm", {
  
 # clustering on individual
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR0")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR0")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR1")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR1")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = "individual", type="CR2")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$state, type = "CR2")[twoway_index,twoway_index])

  # clustering on time
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR0")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR0")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR1")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR1")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = "time", type="CR2")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$year, type = "CR2")[twoway_index,twoway_index])

  # clustering on a randomly generated factor
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR0")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR0")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR1")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR1")[twoway_index,twoway_index])
  expect_equal(vcovCR(plm_twoways, cluster = Produc_scramble$cluster, type="CR2")[twoway_names,twoway_names], 
               vcovCR(lm_twoways, cluster = Produc$cluster, type = "CR2")[twoway_index,twoway_index])
  
})

test_that("vcovCR options work for CR2", {
  n <- nrow(Produc_scramble)
  
  CR2_iv <- vcovCR(plm_individual, type = "CR2")
  expect_identical(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR2"), CR2_iv)
  expect_identical(vcovCR(plm_individual, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_identical(vcovCR(plm_individual, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_individual, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_identical(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_identical(vcovCR(plm_individual, type = "CR2", target = rep(1, n)), CR2_not)
  expect_identical(vcovCR(plm_individual, type = "CR2", target = rep(1, n), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(plm_individual, type = "CR2", target = 1 / Produc_scramble$emp), CR2_not))
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_individual, type = "CR4")
  expect_identical(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR4"), CR4_iv)
  expect_identical(vcovCR(plm_individual, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_identical(vcovCR(plm_individual, type = "CR4", target = rep(1, n), inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_individual, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_identical(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_identical(vcovCR(plm_individual, type = "CR4", target = rep(1, n)), CR4_not)
  expect_identical(vcovCR(plm_individual, type = "CR4", target = rep(1, n), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(plm_individual, type = "CR4", target = 1 / Produc_scramble$emp), CR4_not))
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
  CR_individual <- lapply(CR_types[1:4], function(t) as.matrix(vcovCR(plm_individual, cluster = 1:n, type = t)))
  HC_individual <- lapply(paste0("HC",0:3), function(t) vcovHC(lm_individual, type = t)[individual_index,individual_index])
  expect_equal(CR_individual, HC_individual)
  
  CR_time <- lapply(CR_types[1:4], function(t) as.matrix(vcovCR(plm_time, cluster = 1:n, type = t)))
  HC_time <- lapply(paste0("HC",0:3), function(t) vcovHC(lm_time, type = t)[time_index,time_index])
  expect_equal(CR_time, HC_time)
  
  CR_twoways <- lapply(CR_types[1:4], function(t) as.matrix(vcovCR(plm_twoways, cluster = 1:n, type = t)))
  HC_twoways <- lapply(paste0("HC",0:3), function(t) vcovHC(lm_twoways, type = t)[twoway_index,twoway_index])
  expect_equal(CR_twoways, HC_twoways)
})
