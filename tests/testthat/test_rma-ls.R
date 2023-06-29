context("rma.uni location-scale models")

skip_if_not_installed("metadat")
skip_if_not_installed("metafor")

library(metadat)
suppressMessages(library(metafor, quietly=TRUE))

dat <- dat.bangertdrowns2004
dat$ni100 <- dat$ni/100
dat$meta[is.na(dat$meta)] <- 0
res <- rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data = dat)


test_that("bread works", {
  expect_true(check_bread(res, cluster = dat$id, y = dat$yi))
  vcov_mat <- bread(res) / nobs(res)
  attr(vcov_mat, "dimnames") <- attr(vcov(res)$beta, "dimnames")
  expect_equal(vcov(res)$beta, vcov_mat)
})

CR_types <- paste0("CR",0:4)

test_that("order doesn't matter", {
  
  skip_on_cran()
  
  check_sort_order(res, dat, cluster = "id")
  
})

test_that("clubSandwich works with dropped covariates", {

  dat_miss <- dat.bangertdrowns2004
  expect_warning(res_drop <- rma(yi, vi, mods = ~ length + feedback + info, scale = ~ wic, data = dat))
  
  subset_ind <-  with(dat_miss, complete.cases(length, feedback, info, wic, yi, vi))
  res_complete <- rma(yi, vi, mods = ~ length + feedback + info, scale = ~ wic, data = dat_miss[subset_ind,])
  expect_error(vcovCR(res_complete, type = "CR0", cluster = dat_miss$id))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(res_drop, type = x, cluster = dat_miss$id))
  CR_complete <- lapply(CR_types, function(x) vcovCR(res_complete, type = x, cluster = dat_miss$id[subset_ind]))
  expect_equal(CR_drop, CR_complete)

  test_drop <- lapply(CR_types, function(x) coef_test(res_drop, vcov = x, cluster = dat_miss$id, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(res_complete, vcov = x, cluster = dat_miss$id[subset_ind], test = "All", p_values = FALSE))
  compare_ttests(test_drop, test_complete)
})

test_that("clubSandwich works with missing variances", {
  
  dat_miss <- dat
  dat_miss$vi[sample.int(nrow(dat_miss), size = round(nrow(dat_miss) / 10))] <- NA
  expect_warning(res_drop <- rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data = dat_miss))
  
  
  subset_ind <- with(dat_miss, !is.na(vi))
  res_complete <- rma(yi, vi, mods = ~ ni100 + meta, scale = ~ ni100 + imag, data = dat_miss, subset = !is.na(vi)) 
                       
  expect_error(vcovCR(res_complete, type = "CR0", cluster = dat_miss$id))
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(res_drop, type = x, cluster = dat_miss$id))
  CR_complete <- lapply(CR_types, function(x) vcovCR(res_complete, type = x, cluster = dat_miss$id[subset_ind]))
  expect_equal(CR_drop, CR_complete)
  
})


test_that("vcovCR options work for CR2", {
  RE_var <- res$tau2 + dat$vi
  CR2_iv <- vcovCR(res, type = "CR2", cluster = dat$id)
  expect_equal(vcovCR(res, type = "CR2", cluster = dat$id, inverse_var = TRUE), CR2_iv)

  CR2_not <- vcovCR(res, type = "CR2", cluster = dat$id, inverse_var = FALSE)
  attr(CR2_iv, "inverse_var") <- FALSE
  attr(CR2_iv, "target") <- attr(CR2_not, "target")
  expect_equal(CR2_not, CR2_iv)
  expect_equal(vcovCR(res, type = "CR2", cluster = dat$id, target = RE_var), CR2_not)
  expect_equal(vcovCR(res, type = "CR2", cluster = dat$id, target = RE_var, inverse_var = FALSE), CR2_not)
  expect_false(identical(vcovCR(res, type = "CR2", cluster = dat$id, target = dat$vi), CR2_not))
})
