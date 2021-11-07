context("plm objects - fixed effects")
set.seed(20190513)

library(plm, quietly=TRUE)

data("Produc", package = "plm")
Produc$cluster <- sample(LETTERS[1:10], size = nrow(Produc), replace=TRUE)
Produc_scramble <- Produc[sample(nrow(Produc)),]
n <- nrow(Produc_scramble)

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

test_that("bread works", {
  y <- plm_individual$model$"log(gsp)"
  expect_true(check_bread(plm_individual, cluster = findCluster.plm(plm_individual), y = y))
  sigma_sq_ind <- with(plm_individual, sum(residuals^2) / df.residual) 
  expect_equal(vcov(plm_individual), bread(plm_individual) * sigma_sq_ind / v_scale(plm_individual))
  
  expect_true(check_bread(plm_time, cluster = findCluster.plm(plm_time), y = y))
  sigma_sq_time <- with(plm_time, sum(residuals^2) / df.residual) 
  expect_equal(vcov(plm_time), bread(plm_time) * sigma_sq_time / v_scale(plm_time))
  
  expect_true(check_bread(plm_twoways, cluster = Produc_scramble$state, y = y))
  expect_true(check_bread(plm_twoways, cluster = Produc_scramble$year, y = y))
  sigma_sq_two <- with(plm_twoways, sum(residuals^2) / df.residual) 
  expect_equal(vcov(plm_twoways), bread(plm_twoways) * sigma_sq_two / v_scale(plm_twoways))
})

test_that("CR0 and CR1S agree with arellano vcov", {
  expect_equal(vcovHC(plm_individual, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR0")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_individual, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_individual, type = "CR1S")), check.attributes = FALSE)
  
  expect_equal(vcovHC(plm_time, method="arellano", type = "HC0", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR0")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_time, method="arellano", type = "sss", cluster = "time"), 
               as.matrix(vcovCR(plm_time, type = "CR1S")), check.attributes = FALSE)
  
  expect_equal(vcovHC(plm_twoways, method="arellano", type = "HC0", cluster = "group"), 
               as.matrix(vcovCR(plm_twoways, cluster = "individual", type = "CR0")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_twoways, method="arellano", type = "sss", cluster = "group"), 
               as.matrix(vcovCR(plm_twoways, cluster = "individual", type = "CR1S")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_twoways, method="arellano", type = "HC0", cluster = "time"), 
               as.matrix(vcovCR(plm_twoways, cluster = "time", type = "CR0")), check.attributes = FALSE)
  expect_equal(vcovHC(plm_twoways, method="arellano", type = "sss", cluster = "time"), 
               as.matrix(vcovCR(plm_twoways, cluster = "time", type = "CR1S")), check.attributes = FALSE)
})

test_that("vcovCR options work for CR2", {
  CR2_iv <- vcovCR(plm_individual, type = "CR2")
  expect_equal(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR2"), CR2_iv)
  expect_equal(vcovCR(plm_individual, type = "CR2", inverse_var = TRUE), CR2_iv)
  expect_equal(vcovCR(plm_individual, type = "CR2", target = rep(1, n), inverse_var = TRUE), CR2_iv)
  
  CR2_not <- vcovCR(plm_individual, type = "CR2", inverse_var = FALSE)
  expect_equivalent(CR2_not, CR2_iv)
  expect_equal(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR2", inverse_var = FALSE), CR2_not)
  expect_equal(vcovCR(plm_individual, type = "CR2", target = rep(1, n)), CR2_not)
  expect_equal(vcovCR(plm_individual, type = "CR2", target = rep(1, n), inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(plm_individual, type = "CR2", target = 1 / Produc_scramble$emp), CR2_not))
})

test_that("vcovCR options work for CR4", {
  CR4_iv <- vcovCR(plm_individual, type = "CR4")
  expect_equal(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR4"), CR4_iv)
  expect_equal(vcovCR(plm_individual, type = "CR4", inverse_var = TRUE), CR4_iv)
  expect_equal(vcovCR(plm_individual, type = "CR4", target = rep(1, n), inverse_var = TRUE), CR4_iv)
  
  CR4_not <- vcovCR(plm_individual, type = "CR4", inverse_var = FALSE)
  expect_equivalent(CR4_not, CR4_iv)
  expect_equal(vcovCR(plm_individual, cluster = Produc_scramble$state, type = "CR4", inverse_var = FALSE), CR4_not)
  expect_equal(vcovCR(plm_individual, type = "CR4", target = rep(1, n)), CR4_not)
  expect_equal(vcovCR(plm_individual, type = "CR4", target = rep(1, n), inverse_var = FALSE), CR4_not)
  expect_false(identical(vcovCR(plm_individual, type = "CR4", target = 1 / Produc_scramble$emp), CR4_not))
})

test_that("CR2 and CR4 are target-unbiased", {

  skip_on_cran()
  
  expect_true(check_CR(plm_individual, vcov = "CR2"))
  expect_true(check_CR(plm_individual, vcov = "CR4"))
  expect_true(check_CR(plm_individual, vcov = "CR2", inverse_var = FALSE))
  expect_true(check_CR(plm_individual, vcov = "CR4", inverse_var = FALSE))
  
  expect_true(check_CR(plm_time, vcov = "CR2"))
  expect_true(check_CR(plm_time, vcov = "CR4"))
  expect_true(check_CR(plm_time, vcov = "CR2", inverse_var = FALSE))
  expect_true(check_CR(plm_time, vcov = "CR4", inverse_var = FALSE))
  
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "individual"))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "individual"))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "individual", inverse_var = FALSE))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "individual", inverse_var = FALSE))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "time"))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "time"))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = "time", inverse_var = FALSE))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = "time", inverse_var = FALSE))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = Produc_scramble$cluster))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = Produc_scramble$cluster))
  expect_true(check_CR(plm_twoways, vcov = "CR2", cluster = Produc_scramble$cluster, inverse_var = FALSE))
  expect_true(check_CR(plm_twoways, vcov = "CR4", cluster = Produc_scramble$cluster, inverse_var = FALSE))
})


test_that("vcovCR is equivalent to vcovHC when clusters are all of size 1", {
  library(sandwich, quietly=TRUE)

  CR_types <- paste0("CR",c(0,2))
  HC_types <- paste0("HC",c(0,2))

  CR_individual <- lapply(CR_types, function(t) as.matrix(vcovCR(plm_individual, cluster = 1:n, type = t)))
  HC_individual <- lapply(HC_types, function(t) vcovHC(lm_individual, type = t)[individual_index,individual_index])
  expect_equal(CR_individual, HC_individual)
  
  CR_time <- lapply(CR_types, function(t) as.matrix(vcovCR(plm_time, cluster = 1:n, type = t)))
  HC_time <- lapply(HC_types, function(t) vcovHC(lm_time, type = t)[time_index,time_index])
  expect_equal(CR_time, HC_time)
  
  CR_twoways <- lapply(CR_types, function(t) as.matrix(vcovCR(plm_twoways, cluster = 1:n, type = t)))
  HC_twoways <- lapply(HC_types, function(t) vcovHC(lm_twoways, type = t)[twoway_index,twoway_index])
  expect_equal(CR_twoways, HC_twoways)
})

test_that("CR2 is equivalent to Welch t-test for DiD design", {
  m0 <- 4
  m1 <- 9
  m <- m0 + m1
  cluster <- factor(rep(LETTERS[1:m], each = 2))
  n <- length(cluster)
  time <- rep(c(1,2), m)
  trt_clusters <- c(rep(0,m0), rep(1,m1))
  trt <- (time - 1) * rep(trt_clusters, each = 2)
  nu <- rnorm(m)[cluster]
  e <- rnorm(n)
  y <- 0.4 * trt + nu + e
  
  dat <- data.frame(y, time, trt, cluster)
  plm_DID <- plm(y ~ trt, data = dat, index = c("cluster","time"), 
                 effect = "twoways", model = "within")
  plm_Satt <- coef_test(plm_DID, vcov = "CR2", cluster = dat$cluster)["trt",]
  plm_Wald <- Wald_test(plm_DID, constraints = constrain_zero(1), vcov = "CR2", cluster = dat$cluster)
  df <- m^2 * (m0 - 1) * (m1 - 1) / (m0^2 * (m0 - 1) + m1^2 * (m1 - 1))
  y_diff <- apply(matrix(y, nrow = 2), 2, diff)
  t_Welch <- t.test(y_diff ~ trt_clusters)
  
  expect_equal(with(t_Welch, estimate[[2]] - estimate[[1]]), plm_Satt$beta)
  expect_equal(as.numeric(-t_Welch$statistic), with(plm_Satt, beta / SE))
  expect_equal(as.numeric(-t_Welch$statistic)^2, plm_Wald$Fstat)
  expect_is(all.equal(as.numeric(t_Welch$parameter), plm_Satt$df_Satt), "character")
  expect_equal(plm_Satt$df, df)
  expect_equal(plm_Wald$df_denom, df)
})
