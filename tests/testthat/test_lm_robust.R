context("lm_robust objects")

skip_if_not_installed("estimatr")

library(estimatr)


set.seed(20190513)
data("ChickWeight", package = "datasets")
ChickWeight$wt <- 1 + rpois(nrow(ChickWeight), 3)
ChickWeight$Chick_ordered <- ChickWeight$Chick # James' suggestion 4/16
ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)

lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                    clusters = Chick)
lm_rob_chole <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                          clusters = Chick, try_cholesky = TRUE)

lm_fit_fe <- lm(weight ~ 0 + Time:Diet + Chick, data = ChickWeight)
lm_rob_fe <- lm_robust(weight ~ 0 + Time:Diet, data = ChickWeight, 
                    clusters = Chick, fixed_effects = ~Chick)
lm_rob_fe_chole <- lm_robust(weight ~ 0 + Time:Diet, data = ChickWeight, 
                       clusters = Chick, fixed_effects = ~Chick,
                       try_cholesky = TRUE)

wlm_fit <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, data = ChickWeight)
wlm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                     data = ChickWeight, clusters = Chick)
wlm_rob_chole <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                     data = ChickWeight, clusters = Chick,
                     try_cholesky = TRUE)


# =============== sandwich::bread ===============

test_that("sandwhich::bread works", {
  
  # unweighted tests
  
  bread_lm <- bread(lm_fit)
  bread_rob <- bread(lm_rob)
  bread_lm_fe <- bread(lm_fit_fe)
  bread_rob_fe <- bread(lm_rob_fe)
  
  expect_equal(bread_lm, bread_rob)
  
  # added based on vcovCR tests
  focal_coefs <- names(coef(lm_rob_fe))
  expect_equal(bread_lm_fe[focal_coefs,focal_coefs], as.matrix(bread_rob_fe))
  # expect_equal(bread_lm_fe, bread_rob_fe)
  
  # weighted tests
  
  bread_wlm <- bread(wlm_fit)
  bread_wrob <- bread(wlm_rob)
  # bread_wrob_fe <- bread(wlm_rob_fe)
  
  expect_equal(bread_wlm, bread_wrob)
  # expect_equal(bread_wlm, bread_wrob_fe)
  
})

# =============== model.frame() ===============

test_that("model.frame() works", {
  
  # NOTE: These tests are identical from those for model_matrix(), except, they
  # use model.frame()
  
  # unweighted tests
  
  mf_fit <- model.frame(lm_fit)
  mf_rob <- model.frame(lm_rob)
  mf_rob_chole <- model.frame(lm_rob_chole)
  
  expect_equal(mf_fit, mf_rob)
  expect_equal(mf_fit, mf_rob_chole)
  
  mf_fit_fe <- model.frame(lm_fit_fe)
  mf_rob_fe <- model.frame(lm_rob_fe)
  mf_rob_fe_chole <- model.frame(lm_rob_fe_chole)
  
  expect_equivalent(mf_fit_fe, mf_rob_fe)
  expect_equivalent(mf_fit_fe, mf_rob_fe_chole)
  
  # weighted tests
  
  mf_wlm <- model.frame(wlm_fit)
  mf_wrob <- model.frame(wlm_rob)
  mf_wrob_chole <- model.frame(wlm_rob_chole)
  
  expect_equal(mf_wlm, mf_wrob)
  expect_equal(mf_wlm, mf_wrob_chole)
  
})

# =============== model.matrix() ===============

test_that("model.matrix() works", {
  
  # NOTE: These tests are identical from those for model_matrix(), except, they
  # use model.matrix()
  
  # unweighted tests
  
  mm_fit <- model.matrix(lm_fit)
  mm_rob <- model.matrix(lm_rob)
  mm_rob_chole <- model.matrix(lm_rob_chole)
  
  expect_equal(mm_fit, mm_rob)
  expect_equal(mm_fit, mm_rob_chole)
  
  mm_fit_fe <- model.matrix(lm_fit_fe)
  mm_rob_fe <- augmented_model_matrix(lm_rob_fe)
  mm_rob_fe_chole <- augmented_model_matrix(lm_rob_fe_chole)
  
  expect_equivalent(mm_fit_fe, mm_rob_fe)
  expect_equivalent(mm_fit_fe, mm_rob_fe_chole)
  
  # weighted tests
  
  mm_wlm <- model.matrix(wlm_fit)
  mm_wrob <- model.matrix(wlm_rob)
  mm_wrob_chole <- model.matrix(wlm_rob_chole)
  
  expect_equal(mm_wlm, mm_wrob)
  expect_equal(mm_wlm, mm_wrob_chole)
  
})

# =============== model_matrix() ===============

test_that("model_matrix() works", {
  
  # unweighted tests
  
  mm_fit <- model_matrix(lm_fit) 
  mm_rob <- model_matrix(lm_rob)
  mm_rob_chole <- model_matrix(lm_rob_chole)
  
  expect_equal(mm_fit, mm_rob)
  expect_equal(mm_fit, mm_rob_chole)
  
  mm_fit_fe <- model_matrix(lm_fit_fe)
  mm_rob_fe <- augmented_model_matrix(lm_rob_fe)
  mm_rob_fe_chole <- augmented_model_matrix(lm_rob_fe_chole)
  
  expect_equivalent(mm_fit_fe, mm_rob_fe)
  expect_equivalent(mm_fit_fe, mm_rob_fe_chole)
  
  # weighted tests
  
  mm_wlm <- model_matrix(wlm_fit)
  mm_wrob <- model_matrix(wlm_rob)
  
  expect_equal(mm_wlm, mm_wrob)
  
})

# =============== residuals() ===============

test_that("residuals() works", {
  
  # unweighted tests
  
  res_fit <- residuals(lm_fit)
  res_rob <- residuals(lm_rob)
  res_rob_chole <- residuals(lm_rob_chole)
  
  expect_equal(res_fit, res_rob)
  expect_equal(res_fit, res_rob_chole)
  
  res_fit_fe <- residuals(lm_fit_fe)
  res_rob_fe <- residuals(lm_rob_fe)
  res_rob_fe_chole <- residuals(lm_rob_fe_chole)
  
  expect_equal(res_fit_fe, res_rob_fe)
  expect_equal(res_fit_fe, res_rob_fe_chole)
  
  # weighted tests
  
  res_wlm <- residuals(wlm_fit)
  res_wrob <- residuals(wlm_rob)
  res_wrob_chole <- residuals(wlm_rob_chole)
  
  expect_equal(res_wlm, res_wrob)
  expect_equal(res_wlm, res_wrob_chole)
  
})

# =============== residuals_CS() ===============

test_that("residuals_CS() works", {
  
  # unweighted tests
  
  rcs_fit <- residuals_CS(lm_fit)
  rcs_rob <- residuals_CS(lm_rob)
  rcs_rob_chole <- residuals_CS(lm_rob_chole)
  
  expect_equal(rcs_fit, rcs_rob)
  expect_equal(rcs_fit, rcs_rob_chole)
  
  rcs_fit_fe <- residuals_CS(lm_fit_fe)
  rcs_rob_fe <- residuals_CS(lm_rob_fe)
  rcs_rob_fe_chole <- residuals_CS(lm_rob_fe_chole)
  
  expect_equal(rcs_fit_fe, rcs_rob_fe)
  expect_equal(rcs_fit_fe, rcs_rob_fe_chole)
  
  # weighted tests
  
  rcs_wlm <- residuals_CS(wlm_fit)
  rcs_wrob <- residuals_CS(wlm_rob)
  rcs_wrob_chole <- residuals_CS(wlm_rob_chole)
  
  expect_equal(rcs_wlm, rcs_wrob)
  expect_equal(rcs_wlm, rcs_wrob_chole)
})

# =============== coef() ===============

test_that("coef() works", {
  
  # unweighted tests
  
  coef_fit <- coef(lm_fit)
  coef_rob <- coef(lm_rob)
  coef_rob_chole <- coef(lm_rob_chole)
  
  expect_equal(coef_fit, coef_rob)
  expect_equal(coef_fit, coef_rob_chole)
  
  coef_fit_fe <- coef(lm_fit_fe)
  coef_rob_fe <- coef(lm_rob_fe)
  coef_rob_fe_chole <- coef(lm_rob_fe_chole)
  
  expect_equal(coef_fit_fe[names(coef_rob_fe)], coef_rob_fe)
  expect_equal(coef_fit_fe[names(coef_rob_fe)], coef_rob_fe_chole)
  
  # weighted tests
  
  coef_wlm <- coef(wlm_fit)
  coef_wrob <- coef(wlm_rob)
  coef_wrob_chole <- coef(wlm_rob_chole)
  
  expect_equal(coef_wlm, coef_wrob)
  expect_equal(coef_wlm, coef_wrob_chole)
  
})

# =============== nobs() ===============

test_that("nobs() works", {
  
  # unweighted tests
  
  nobs_fit <- nobs(lm_fit)
  nobs_rob <- nobs(lm_rob)
  nobs_rob_chole <- nobs(lm_rob_chole)
  
  expect_equal(nobs_fit, nobs_rob)
  expect_equal(nobs_fit, nobs_rob_chole)
  
  nobs_rob_fe <- nobs(lm_rob_fe)
  nobs_rob_fe_chole <- nobs(lm_rob_fe_chole)
  
  expect_equal(nobs_fit, nobs_rob_fe)
  expect_equal(nobs_fit, nobs_rob_fe_chole)
  
  # weighted tests
  
  nobs_wlm <- nobs(wlm_fit)
  nobs_wrob <- nobs(wlm_rob)
  nobs_wrob_chole <- nobs(wlm_rob_chole)
  
  expect_equal(nobs_wlm, nobs_wrob)
  expect_equal(nobs_wlm, nobs_wrob_chole)
  
})

# =============== targetVariance() ===============

test_that("targetVariance() works", {
  
  # unweighted tests
  
  tV_fit <- targetVariance(lm_fit, ChickWeight$Chick)
  tV_rob <- targetVariance(lm_rob, ChickWeight$Chick)
  tV_rob_chole <- targetVariance(lm_rob_chole, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob)
  expect_equal(tV_fit, tV_rob_chole)
  
  tV_rob_fe <- targetVariance(lm_rob_fe, ChickWeight$Chick)
  tV_rob_fe_chole <- targetVariance(lm_rob_fe_chole, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob_fe)
  expect_equal(tV_fit, tV_rob_fe_chole)
  
  # weighted tests
  
  tV_wlm <- targetVariance(wlm_fit, ChickWeight$Chick)
  tV_wrob <- targetVariance(wlm_rob, ChickWeight$Chick)
  tV_wrob_chole <- targetVariance(wlm_rob_chole, ChickWeight$Chick)
  
  expect_equal(tV_wlm, tV_wrob)
  expect_equal(tV_wlm, tV_wrob_chole)
  
})

# =============== weightMatrix() ===============

test_that("weightMatrix() works", {
  
  # unweighted tests
  
  wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick)
  wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick)
  wM_rob_chole <- weightMatrix(lm_rob_chole, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob)
  expect_equal(wM_fit, wM_rob_chole)
  
  wM_rob_fe <- weightMatrix(lm_rob_fe, ChickWeight$Chick)
  wM_rob_fe_chole <- weightMatrix(lm_rob_fe_chole, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob_fe)
  expect_equal(wM_fit, wM_rob_fe_chole)
  
  # weighted tests
  
  wM_wlm <- weightMatrix(wlm_fit, ChickWeight$Chick)
  wM_wrob <- weightMatrix(wlm_rob, ChickWeight$Chick)
  wM_wrob_chole <- weightMatrix(wlm_rob_chole, ChickWeight$Chick)
  
  expect_equal(wM_wlm, wM_wrob)
  expect_equal(wM_wlm, wM_wrob_chole)
  
})

# =============== v_scale() ===============

test_that("v_scale() works", {
  
  # unweighted tests
  
  vs_fit <- v_scale(lm_fit)
  vs_rob <- v_scale(lm_rob)
  vs_rob_chole <- v_scale(lm_rob_chole)
  vs_rob_fe <- v_scale(lm_rob_fe)
  vs_rob_fe_chole <- v_scale(lm_rob_fe_chole)
  
  expect_equal(vs_fit, vs_rob)
  expect_equal(vs_fit, vs_rob_chole)
  expect_equal(vs_fit, vs_rob_fe)
  expect_equal(vs_fit, vs_rob_fe_chole)
  
  # weighted tests
  
  vs_wlm <- v_scale(wlm_fit)
  vs_wrob <- v_scale(wlm_rob)
  vs_wrob_chole <- v_scale(wlm_rob_chole)
  
  expect_equal(vs_wlm, vs_wrob)
  expect_equal(vs_wlm, vs_wrob_chole)
  
})

# =============== vcovCR ===============

test_that("vcovCR works", {
  
  types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")
  
  # unweighted tests
  
  focal_coefs <- names(coef(lm_rob_fe))
  
  for (type in types) {
    
    vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, type = type)
    vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, type = type)
    vcov_lmr_chole <- vcovCR(lm_rob_chole, ChickWeight$Chick, type = type)
    
    expect_equal(vcov_lm, vcov_lmr, 
                 label = paste0("When type = ", type, ", ", "vcov_lm"))
    expect_equal(vcov_lm, vcov_lmr_chole, 
                 label = paste0("When type = ", type, ", ", "vcov_lm"))
    
    if (type %in% c("CR0", "CR1", "CR2")) {
      
      vcov_lm_fe <- vcovCR(lm_fit_fe, ChickWeight$Chick, type = type)
      vcov_lmr_fe <- vcovCR(lm_rob_fe, ChickWeight$Chick, type = type)
      vcov_lmr_fe_chole <- vcovCR(lm_rob_fe_chole, ChickWeight$Chick, type = type)
      
      expect_equal(vcov_lm_fe[focal_coefs,focal_coefs], as.matrix(vcov_lmr_fe), 
                   label = paste0("When type = ", type, ", ", "vcov_lm_fe[focal_coefs,focal_coefs]"))
      expect_equal(vcov_lm_fe[focal_coefs,focal_coefs], as.matrix(vcov_lmr_fe_chole), 
                   label = paste0("When type = ", type, ", ", "vcov_lm_fe[focal_coefs,focal_coefs]"))
      
    }
    
    # weighted tests
    
    vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, type = type)
    vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, type = type)
    vcov_wlmr_chole <- vcovCR(wlm_rob_chole, ChickWeight$Chick, type = type)
    
    expect_equal(vcov_wlm, vcov_wlmr)
    expect_equal(vcov_wlm, vcov_wlmr_chole)
    
    
    if (type %in% c("CR0","CR2")) {
      
      lm_rob_type <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, 
        se_type = type
      )
      lm_rob_type_chole <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, 
        se_type = type,
        try_cholesky = TRUE
      )
      
      expect_equal(as.matrix(vcov_lmr), vcov(lm_rob_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr)"))
      expect_equal(as.matrix(vcov_lmr), vcov(lm_rob_type_chole), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr)"))
      
      lm_rob_fe_type <- lm_robust(
        weight ~ 0 + Time:Diet, data = ChickWeight, 
        clusters = Chick, fixed_effects = ~Chick,
        se_type = type
      )
      lm_rob_fe_type_chole <- lm_robust(
        weight ~ 0 + Time:Diet, data = ChickWeight, 
        clusters = Chick, fixed_effects = ~Chick,
        se_type = type,
        try_cholesky = TRUE
      )
      
      expect_equal(as.matrix(vcov_lmr_fe), vcov(lm_rob_fe_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr_fe)"))
      expect_equal(as.matrix(vcov_lmr_fe), vcov(lm_rob_fe_type_chole), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr_fe)"))
      
      wlm_rob_type <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, weights = wt,
        se_type = type
      )
      wlm_rob_type_chole <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, weights = wt,
        se_type = type,
        try_cholesky = TRUE
      )
      
      expect_equal(as.matrix(vcov_wlmr), vcov(wlm_rob_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_wlmr)"))
      expect_equal(as.matrix(vcov_wlmr), vcov(wlm_rob_type_chole), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_wlmr)"))
      
    }
  }
  
})

test_that("vovCR properly pulls cluster specified for lm_robust", {
  
  # unweighted tests
  
  uw_clust <- vcovCR(lm_rob, ChickWeight$Chick, "CR2")
  uw_no_clust <- vcovCR(lm_rob, type = "CR2")
  uw_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2")
  
  expect_equal(uw_clust, uw_no_clust)
  expect_equal(uw_no_clust, uw_lm)
  
  # create an lm_robust that draws in data differently
  lm_rob_fact <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                           clusters = factor(ChickWeight$Chick_ordered, ordered = FALSE))
  # perform vcovCR
  uw_fact_cr <- vcovCR(lm_rob_fact, type = "CR2")
  
  # check they are the same
  expect_equivalent(uw_clust, uw_fact_cr)
  
  # put cluster data in a variable
  # fact <- factor(ChickWeight$Chick_ordered, ordered = FALSE)
  fact <- ChickWeight$Chick
  
  # pass variable to lm_robust
  lm_rob_var <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                           clusters = fact)
  
  # perform vcovCR
  uw_fact_var <- vcovCR(lm_rob_var, type = "CR2")
  
  # check they are the same
  expect_equivalent(uw_clust, uw_fact_var)
  
  # weighted tests
  
  w_clust <- vcovCR(wlm_rob, ChickWeight$Chick, "CR2")
  w_no_clust <- vcovCR(wlm_rob, type = "CR2")
  w_lm <- vcovCR(wlm_fit, ChickWeight$Chick, "CR2")
  
  expect_equal(w_clust, w_no_clust)
  expect_equal(w_no_clust, w_lm)
  
  # create an lm_robust that draws in data differently
  lm_rob_fact_w <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                             data = ChickWeight, 
                             clusters = factor(ChickWeight$Chick_ordered, ordered = FALSE))
  # perform vcovCR
  w_fact_cr <- vcovCR(lm_rob_fact_w, type = "CR2")
  
  expect_equal(w_clust, w_fact_cr)
  
  # pass variable to lm_robust
  lm_rob_var_w <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt,
                            data = ChickWeight, clusters = fact)
  
  # perform vcovCR
  w_fact_var <- vcovCR(lm_rob_var_w, type = "CR2")
  
  # check they are the same
  expect_equal(w_clust, w_fact_var)
  
})

# =============== na.action.lm_robust() ===============

test_that("na.action.lm_robust() works correctly", {
  
  compare_na_actions <- function(i) {
    
    # generate random data
    n <- 100
    df <- data.frame(
      y = rnorm(n),
      x1 = rnorm(n),
      x2 = rnorm(n)
    )
    
    # add random NA values to y and x2
    miss_rows <- sample.int(n, size = n/10)
    df$y[miss_rows] <- NA
    df$x2[miss_rows] <- NA
    
    # fit models
    linear <- lm(y ~ x1 + x2 + x1:x2, data = df)
    robust <- lm_robust(y ~ x1 + x2 + x1:x2, data = df)
    
    # get na.action() of models
    na_lm <- na.action(linear)
    na_rob <- na.action(robust)
    
    # compare
    expect_equal(na_lm, na_rob)
  }
  
  # compare 10 times with different random data
  lapply(1:10, compare_na_actions)
  
})

# =============== Higher level Tests ===============


data("Seatbelts", package = "datasets")

# Convert Seatbelts time series to data frame
belts <- as.data.frame(Seatbelts)

# Extract the time index and convert to Date
time_index <- time(Seatbelts)
year <- floor(time_index)
month <- cycle(Seatbelts)

# Add the time columns
belts$kms <- belts$kms - mean(belts$kms)
belts$year <- year - mean(year)
belts$month <- month

# Create identical lm and lm_robust models
belts_fit <- lm(DriversKilled ~ kms + PetrolPrice + law + year, data = belts)
belts_rob <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts)


test_that("test Wald_test() with lm_robust", {
  
  Wald_FIT <- Wald_test(
    belts_fit,
    constraints = constrain_zero("year", reg_ex = TRUE),
    vcov = "CR2",
    cluster = belts$month
  )
  
  Wald_ROB <- Wald_test(
    belts_rob,
    constraints = constrain_zero("year", reg_ex = TRUE),
    vcov = "CR2",
    cluster = belts$month
  )
  
  expect_equal(Wald_FIT, Wald_ROB)
  
})


test_that("test conf_int() with lm_robust", {
  
  conf_FIT <- conf_int(belts_fit, vcov = "CR2", cluster = belts$month)
  conf_ROB <- conf_int(belts_rob, vcov = "CR2", cluster = belts$month)
  
  expect_equal(conf_FIT, conf_ROB)
  
})


test_that("test coef_test() with lm_robust", {
  
  coef_FIT <- coef_test(belts_fit, vcov = "CR2", cluster = belts$month)
  coef_ROB <- coef_test(belts_rob, vcov = "CR2", cluster = belts$month)
  
  expect_equal(coef_FIT, coef_ROB)
  
})


# =============== Tests Based on test_lm.R ===============


test_that("Order doesn't matter.",{
  
  check_sort_order(belts_rob, belts, "month")
  
})


test_that("clubSandwich works with dropped observations", {
  belts_miss <- belts
  miss_indicator <- sample.int(nrow(belts), size = round(nrow(belts) / 10))
  belts_miss$law[miss_indicator] <- NA
  belts_miss$kms[miss_indicator] <- NA
  
  rob_dropped <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts_miss, cluster = month)
  belts_complete <- subset(belts_miss, !is.na(law))
  rob_complete <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts_complete, cluster = month)
  
  CR_types <- paste0("CR",0:4)
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(rob_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(rob_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(rob_dropped, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(rob_complete, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


test_that("clubSandwich requires no missing values on the clustering variable", {
  belts_miss <- belts
  miss_indicator <- sample.int(nrow(belts), size = round(nrow(belts) / 10))
  belts_miss$month[miss_indicator] <- NA
  
  rob_dropped <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts_miss)
  
  expect_error(vcovCR(rob_dropped, cluster = belts_miss$month, type = "CR0"), 
               "Clustering variable cannot have missing values.")
  expect_error(coef_test(rob_dropped, vcov = "CR0", cluster = belts_miss$month, test = "All"),
               "Clustering variable cannot have missing values.")
})
