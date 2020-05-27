context("multi-variate multi-level lme objects")
skip()
set.seed(20190513)

dat <- read.table(file="https://raw.githubusercontent.com/wviechtb/multivariate_multilevel_models/master/data.dat", header=TRUE, sep="\t")
dat$pa <- rowMeans(dat[, grepl("pa", names(dat))])
dat$na <- rowMeans(dat[, grepl("na", names(dat))])

# keep only variables that are needed
dat <- dat[, c("id", "sex", "beep", "pa", "na")]

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

obj <- MVML_full
cluster <- dat$id
y <- dat$y
targetVariance(obj)

test_that("bread works", {
  bread_checks <- lapply(objects, check_bread, cluster = dat$id, y = dat$y)
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
  
})


test_that("CR2 is target-unbiased", {
  skip_on_cran()
  CR2_checks <- mapply(check_CR, obj = objects, vcov = CR2_mats)
  expect_true(all(CR2_checks))
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

  obj <- obj_dropped
  cluster <- nlme::getGroups(obj, level = 1)
  target <- NULL
  inverse_var <- is.null(target)
  type <- "CR2"
  form <- "sandwich"
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(obj_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(obj_complete, type = x))
  expect_identical(CR_drop, CR_complete)

  test_drop <- lapply(CR_drop, function(x) coef_test(obj_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_complete, function(x) coef_test(obj_complete, vcov = x, test = "All", p_values = FALSE))
  expect_identical(test_drop, test_complete)
})
