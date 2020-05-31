context("3-level lme objects")
set.seed(20190513)

suppressMessages(library(lme4, quietly=TRUE))
library(nlme, quietly=TRUE, warn.conflicts=FALSE)
library(mlmRev, quietly=TRUE, warn.conflicts=FALSE)

school_subset <- levels(egsingle$schoolid)
school_subset <- sample(school_subset, size = 15)
egsingle <- droplevels(subset(egsingle, schoolid %in% school_subset))


obj_A1 <- lme(math ~ year * size + female + black + hispanic,
              random = list(~ year | schoolid, ~ 1 | childid),
              data = egsingle)
obj_A2 <- update(obj_A1, weights = varIdent(form = ~ 1 | female))
obj_A3 <- update(obj_A1, correlation = corExp(form = ~ year))
obj_A4 <- update(obj_A2, correlation = corExp(form = ~ year))
objects <- list(A1 = obj_A1, A2 = obj_A2, A3 = obj_A3, A4 = obj_A4)

CR2_mats <- lapply(objects, vcovCR, type = "CR2")

test_that("bread works", {
  bread_checks <- lapply(objects, check_bread, cluster = egsingle$schoolid, y = egsingle$math)
  expect_true(all(unlist(bread_checks)))
  
  obj_vcovs <- lapply(objects, vcov)
  obj_bread <- lapply(objects, function(obj) obj$sigma^2 * sandwich::bread(obj) / v_scale(obj))
  expect_equal(obj_vcovs, obj_bread)
})


test_that("vcovCR options work for CR2", {
  skip_on_cran()
  
  expect_identical(vcovCR(obj_A1, cluster = egsingle$schoolid, type = "CR2"), CR2_mats[["A1"]])
  expect_equal(vcovCR(obj_A1, type = "CR2", inverse_var = TRUE), CR2_mats[["A1"]])
  expect_false(identical(vcovCR(obj_A1, type = "CR2", inverse_var = FALSE), CR2_mats[["A1"]]))
  target <- targetVariance(obj_A1)
  expect_equal(vcovCR(obj_A1, type = "CR2", target = target, inverse_var = TRUE), CR2_mats[["A1"]])
  attr(CR2_mats[["A1"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A1, type = "CR2", target = target, inverse_var = FALSE), CR2_mats[["A1"]])

  expect_identical(vcovCR(obj_A2, cluster = egsingle$schoolid, type = "CR2"), CR2_mats[["A2"]])
  expect_equal(vcovCR(obj_A2, type = "CR2", inverse_var = TRUE), CR2_mats[["A2"]])
  expect_false(identical(vcovCR(obj_A2, type = "CR2", inverse_var = FALSE), CR2_mats[["A2"]]))
  target <- targetVariance(obj_A2)
  expect_equal(vcovCR(obj_A2, type = "CR2", target = target, inverse_var = TRUE), CR2_mats[["A2"]])
  attr(CR2_mats[["A2"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A2, type = "CR2", target = target, inverse_var = FALSE), CR2_mats[["A2"]])

  expect_identical(vcovCR(obj_A3, cluster = egsingle$schoolid, type = "CR2"), CR2_mats[["A3"]])
  expect_equal(vcovCR(obj_A3, type = "CR2", inverse_var = TRUE), CR2_mats[["A3"]])
  expect_false(identical(vcovCR(obj_A3, type = "CR2", inverse_var = FALSE), CR2_mats[["A3"]]))
  target <- targetVariance(obj_A3)
  expect_equal(vcovCR(obj_A3, type = "CR2", target = target, inverse_var = TRUE), CR2_mats[["A3"]])
  attr(CR2_mats[["A3"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A3, type = "CR2", target = target, inverse_var = FALSE), CR2_mats[["A3"]])
  
  expect_identical(vcovCR(obj_A4, cluster = egsingle$schoolid, type = "CR2"), CR2_mats[["A4"]])
  expect_equal(vcovCR(obj_A4, type = "CR2", inverse_var = TRUE), CR2_mats[["A4"]])
  expect_false(identical(vcovCR(obj_A4, type = "CR2", inverse_var = FALSE), CR2_mats[["A4"]]))
  target <- targetVariance(obj_A4)
  expect_equal(vcovCR(obj_A4, type = "CR2", target = target, inverse_var = TRUE), CR2_mats[["A4"]])
  attr(CR2_mats[["A4"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A4, type = "CR2", target = target, inverse_var = FALSE), CR2_mats[["A4"]])
  
})

test_that("vcovCR options work for CR4", {
  skip_on_cran()
  skip("Not worrying about CR4 for now.")
  CR4_mats <- lapply(objects, vcovCR, type = "CR4")
  
  expect_identical(vcovCR(obj_A1, cluster = egsingle$schoolid, type = "CR4"), CR4_mats[["A1"]])
  expect_identical(vcovCR(obj_A1, type = "CR4", inverse_var = TRUE), CR4_mats[["A1"]])
  expect_false(identical(vcovCR(obj_A1, type = "CR4", inverse_var = FALSE), CR4_mats[["A1"]]))
  target <- targetVariance(obj_A1)
  expect_equal(vcovCR(obj_A1, type = "CR4", target = target, inverse_var = TRUE), CR4_mats[["A1"]])
  attr(CR4_mats[["A1"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A1, type = "CR4", target = target, inverse_var = FALSE), CR4_mats[["A1"]])
  
  expect_identical(vcovCR(obj_A2, cluster = egsingle$schoolid, type = "CR4"), CR4_mats[["A2"]])
  expect_identical(vcovCR(obj_A2, type = "CR4", inverse_var = TRUE), CR4_mats[["A2"]])
  expect_false(identical(vcovCR(obj_A2, type = "CR4", inverse_var = FALSE), CR4_mats[["A2"]]))
  target <- targetVariance(obj_A2)
  expect_equal(vcovCR(obj_A2, type = "CR4", target = target, inverse_var = TRUE), CR4_mats[["A2"]])
  attr(CR4_mats[["A2"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A2, type = "CR4", target = target, inverse_var = FALSE), CR4_mats[["A2"]])
  
  expect_identical(vcovCR(obj_A3, cluster = egsingle$schoolid, type = "CR4"), CR4_mats[["A3"]])
  expect_identical(vcovCR(obj_A3, type = "CR4", inverse_var = TRUE), CR4_mats[["A3"]])
  expect_false(identical(vcovCR(obj_A3, type = "CR4", inverse_var = FALSE), CR4_mats[["A3"]]))
  target <- targetVariance(obj_A3)
  expect_equal(vcovCR(obj_A3, type = "CR4", target = target, inverse_var = TRUE), CR4_mats[["A3"]])
  attr(CR4_mats[["A3"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A3, type = "CR4", target = target, inverse_var = FALSE), CR4_mats[["A3"]])
  
  expect_identical(vcovCR(obj_A4, cluster = egsingle$schoolid, type = "CR4"), CR4_mats[["A4"]])
  expect_identical(vcovCR(obj_A4, type = "CR4", inverse_var = TRUE), CR4_mats[["A4"]])
  expect_false(identical(vcovCR(obj_A4, type = "CR4", inverse_var = FALSE), CR4_mats[["A4"]]))
  target <- targetVariance(obj_A4)
  expect_equal(vcovCR(obj_A4, type = "CR4", target = target, inverse_var = TRUE), CR4_mats[["A4"]])
  attr(CR4_mats[["A4"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(obj_A4, type = "CR4", target = target, inverse_var = FALSE), CR4_mats[["A4"]])
  
})


test_that("CR2 is target-unbiased", {
  skip_on_cran()
  CR2_checks <- mapply(check_CR, obj = objects, vcov = CR2_mats)
  expect_true(all(CR2_checks))
  # CR4_checks <- mapply(check_CR, obj = objects, vcov = CR4_mats)
  # expect_true(all(CR4_checks))
})


CR_types <- paste0("CR",0:3)

test_that("Order doesn't matter.", {
  skip_on_cran()
  
  check_sort_order(obj_A4, egsingle)
  
})


test_that("clubSandwich works with dropped observations", {
  skip_on_cran()
  dat_miss <- egsingle
  dat_miss$math[sample.int(nrow(egsingle), size = round(nrow(egsingle) / 10))] <- NA
  obj_dropped <- update(obj_A4, data = dat_miss, na.action = na.omit)
  obj_complete <- update(obj_A4, data = dat_miss, subset = !is.na(math))

  # obj <- obj_dropped
  # cluster <- nlme::getGroups(obj, level = 1)
  # target <- NULL
  # inverse_var <- TRUE
  # type <- "CR2"
  # form <- "sandwich"
  # 
  # full_grps <- get_cor_grouping(obj)
  # R_list <- nlme::corMatrix(obj$modelStruct$corStruct)
  # levels <- names(R_list)
  # grps <- get_cor_grouping(obj, levels = names(R_list))
  # 
  # V_list <- build_var_cor_mats(obj)
  # V_grps <- attr(V_list, "groups")
  # ZDZ_list <- build_RE_mats(obj)
  # ZDZ_grps <- attr(ZDZ_list, "groups")
  # 
  # V_dim <- sapply(V_list, nrow)
  # identical(names(V_dim), names(table(V_grps)))
  # data.frame(dim = V_dim, grps = table(V_grps))
  # table(V_dim == table(V_grps))
  # dat_miss$x <- NA
  # dat_miss$x[!is.na(dat_miss$math)] <- V_grps
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_identical(CR_drop, CR_complete)

  test_drop <- lapply(CR_drop, function(x) coef_test(obj_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_complete, function(x) coef_test(obj_complete, vcov = x, test = "All", p_values = FALSE))
  expect_identical(test_drop, test_complete)
})


test_that("Possible to cluster at higher level than random effects", {
  skip_on_cran()
  
  # fit two-level model
  obj_2level <- lme(math ~ year * size + female + black + hispanic,
                    random = ~ year | childid,
                    data = egsingle)
  
  # cluster at level 3
  V <- vcovCR(obj_2level, type = "CR2", cluster = egsingle$schoolid)
  expect_is(V, "vcovCR")
  
  # create 4th level
  n_districts <- nlevels(egsingle$schoolid) / 3
  districtid <- rep(1:n_districts, each = 3)[egsingle$schoolid]
  
  # cluster at level 4
  expect_is(vcovCR(obj_2level, type = "CR2", cluster = districtid), "vcovCR")
  expect_is(vcovCR(obj_A1, type = "CR2", cluster = districtid), "vcovCR")
  expect_is(vcovCR(obj_A2, type = "CR2", cluster = districtid), "vcovCR")
  expect_is(vcovCR(obj_A3, type = "CR2", cluster = districtid), "vcovCR")
  expect_is(vcovCR(obj_A4, type = "CR2", cluster = districtid), "vcovCR")
})
