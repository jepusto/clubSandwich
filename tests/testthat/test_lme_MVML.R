context("multi-variate multi-level lme objects")
set.seed(20190513)

skip_on_cran()

dat <- read.table(file="https://raw.githubusercontent.com/wviechtb/multivariate_multilevel_models/master/data.dat", header=TRUE, sep="\t")
dat$pa <- rowMeans(dat[, grepl("pa", names(dat))])
dat$na <- rowMeans(dat[, grepl("na", names(dat))])

# keep only variables that are needed
dat <- dat[, c("id", "sex", "beep", "pa", "na")]

# keep only the first 10 IDs
dat <- subset(dat, id <= 10)

# change into very long format
dat <- reshape(dat, direction="long", varying=c("pa","na"), v.names="y", idvar="obs", timevar="outcome")
dat$obs <- NULL
dat <- dat[order(dat$id, dat$beep, dat$outcome),]
rownames(dat) <- 1:nrow(dat)
dat$outnum <- dat$outcome
dat$outcome <- factor(dat$outcome)


library(nlme, quietly=TRUE, warn.conflicts=FALSE)

MVML_full <- lme(y ~ outcome - 1,
                 random = ~ outcome - 1 | id,
                 weights = varIdent(form = ~ 1 | outcome),
                 correlation = corSymm(form = ~ outnum | id/beep),
                 data = dat, na.action = na.omit)

MVML_diag <- lme(y ~ outcome - 1,
                 random = ~ outcome - 1 | id,
                 weights = varIdent(form = ~ 1 | outcome),
                 data = dat, na.action = na.omit)

gls_full <- gls(y ~ outcome - 1,
                weights = varIdent(form = ~ 1 | outcome),
                correlation = corSymm(form = ~ outnum | id/beep),
                data = dat, na.action = na.omit)

objects <- list(MVML_full = MVML_full, MVML_diag = MVML_diag, gls = gls_full)

CR2_mats <- lapply(objects, vcovCR, type = "CR2")

test_that("bread works", {
  bread_checks <- lapply(objects, check_bread, cluster = dat$id, y = dat$y)
  expect_true(all(unlist(bread_checks)))
  
  obj_vcovs <- lapply(objects, vcov)
  obj_bread <- lapply(objects, function(obj) obj$sigma^2 * sandwich::bread(obj) / v_scale(obj))
  expect_equal(obj_vcovs, obj_bread)
})


test_that("vcovCR options work for CR2", {
  
  expect_identical(vcovCR(MVML_full, cluster = dat$id, type = "CR2"), CR2_mats[["MVML_full"]])
  expect_equal(vcovCR(MVML_full, type = "CR2", inverse_var = TRUE), CR2_mats[["MVML_full"]])
  expect_false(identical(vcovCR(MVML_full, type = "CR2", inverse_var = FALSE), CR2_mats[["MVML_full"]]))
  target <- targetVariance(MVML_full)
  expect_equal(vcovCR(MVML_full, type = "CR2", target = target, inverse_var = TRUE), CR2_mats[["MVML_full"]])
  attr(CR2_mats[["MVML_full"]], "inverse_var") <- FALSE
  expect_equal(vcovCR(MVML_full, type = "CR2", target = target, inverse_var = FALSE), CR2_mats[["MVML_full"]])
  
})


test_that("CR2 is target-unbiased", {
  
  CR2_checks <- mapply(check_CR, obj = objects, vcov = CR2_mats)
  expect_true(all(CR2_checks))
})


CR_types <- paste0("CR",0:3)

test_that("Order doesn't matter.", {
  
  check_sort_order(MVML_full, dat, seed = 20200530)
  check_sort_order(MVML_diag, dat, seed = 20200530)
  
})


test_that("clubSandwich works with dropped observations", {
  
  dat_miss <- dat
  dat_miss$y[sample.int(nrow(dat), size = round(nrow(dat) / 10))] <- NA
  obj_dropped <- update(MVML_full, data = dat_miss, na.action = na.omit)
  obj_complete <- update(MVML_full, data = dat_miss, subset = !is.na(y))

  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_identical(CR_drop, CR_complete)

  test_drop <- lapply(CR_drop, function(x) coef_test(obj_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_complete, function(x) coef_test(obj_complete, vcov = x, test = "All", p_values = FALSE))
  expect_identical(test_drop, test_complete)
})


test_that("Possible to cluster at higher level than random effects", {
  
  # create 4th level
  n_groups <- nlevels(factor(dat$id))
  group_id <- rep(1:n_groups, each = 4)[dat$id]
  
  # cluster at level 4
  expect_is(vcovCR(MVML_full, type = "CR2", cluster = group_id), "vcovCR")
  expect_is(vcovCR(MVML_diag, type = "CR2", cluster = group_id), "vcovCR")
  
})
