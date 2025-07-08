context("lm_lin robust objects")

skip_if_not_installed("estimatr")

suppressMessages(library(estimatr, quietly=TRUE))


set.seed(20250630)
data("Baumann", package = "carData")
Baumann$pretest.1_c <- Baumann$pretest.1 - mean(Baumann$pretest.1)
Baumann$pretest.2_c <- Baumann$pretest.2 - mean(Baumann$pretest.2)

lin <- lm_lin(
  post.test.3 ~ group,
  covariates = ~ pretest.1 + pretest.2,
  data = Baumann
)
rob <- lm_robust(
  post.test.3 ~ group + pretest.1_c + pretest.2_c + 
    group:pretest.1_c + group:pretest.2_c,
  data = Baumann
)


data("ChickWeight", package = "datasets")
N <- nrow(ChickWeight)
ChickWeight$wt <- 1 + rpois(N, 3)
ChickWeight$Chick_ordered <- ChickWeight$Chick # James' suggestion 4/16
ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)
ChickWeight$Chick_int <- as.integer(ChickWeight$Chick)
ChickWeight$rando <- "Drop"
ChickWeight$rando[sample(1:N, size = round(0.8 * N))] <- "Keep"
table(ChickWeight$rando)

lin_chick <- lm_lin(
  weight ~ Diet,
  covariates = ~ Time,
  data = ChickWeight,
  clusters = Chick
)

rob_chick <- lm_robust(
  weight ~ Diet + Time + Diet:Time,
  data = ChickWeight,
  clusters = Chick
)

test_that("model.frame() works", {
  
  mf_lin <- model.frame(lin)
  mf_rob <- model.frame(rob)
  
  # focal_cols <- names(mf_lin)
  
  expect_equivalent(mf_lin, mf_rob)
  
})

test_that("model_matrix() works", {
  
  mm_lin <- model_matrix(lin)
  mm_rob <- model_matrix(rob)
  
  expect_equivalent(mm_lin, mm_rob)
  
})

test_that("model_matrix() works without fixest", {
  
  local_mocked_bindings(
    requireNamespace = function(...) FALSE
  )
  
  expect_false(requireNamespace("fixest", quietly = TRUE))
  
  # unweighted tests
  
  mm_lin <- model_matrix(lin)
  mm_rob <- model_matrix(rob)
  # amm_rob_fe <- augmented_model_matrix(rob_fe)
  
  expect_equivalent(mm_lin, mm_rob)
  
  # mm_fit_fe <- model_matrix(lin_fe)
  # mm_rob_fe <- model_matrix(rob_fe)

  # Check that fixed effects are the same
  # expect_equivalent(mm_fit_fe[,colnames(amm_rob_fe)], amm_rob_fe)
  # Core predictor matrices are different
  # expect_false(identical(mm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  # But one can be computed by residualizing
  # expect_equivalent(
  #   mm_rob_fe,
  #   residuals(lm.fit(amm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  # )
  # model matrix is centered by chick so means are zero
  # expect_lt(max(abs(apply(mm_rob_fe, 2, \(x) tapply(x, ChickWeight$Chick, mean)))), 1e-12)
  
  # expect_message(
    # vcovCR(rob_fe,type = "CR0"),
    # "For improved performance in models with fixed effects, install the package \\{fixest\\}\\."
  # )
  
})

test_that("targetVariance() works", {

  # unweighted tests

  tV_fit <- targetVariance(lin_chick, ChickWeight$Chick)
  tV_rob <- targetVariance(rob_chick, ChickWeight$Chick)

  expect_equal(tV_fit, tV_rob)
#   
#   tV_rob_fe <- targetVariance(rob_fe, ChickWeight$Chick)
#   
#   expect_equal(tV_fit, tV_rob_fe)
#   
#   # weighted tests
#   
#   tV_wlm <- targetVariance(wlin, ChickWeight$Chick)
#   tV_wrob <- targetVariance(wrob, ChickWeight$Chick)
#   
#   expect_equal(tV_wlm, tV_wrob)
#   
})


test_that("weightMatrix() works", {

  # unweighted tests

  wM_fit <- weightMatrix(lin_chick, ChickWeight$Chick)
  wM_rob <- weightMatrix(rob_chick, ChickWeight$Chick)

  expect_equal(wM_fit, wM_rob)

  # weighted tests
#   
#   wM_wlm <- weightMatrix(wlin, ChickWeight$Chick)
#   wM_wrob <- weightMatrix(wrob, ChickWeight$Chick)
#   
#   expect_equal(wM_wlm, wM_wrob)
#   
})

test_that("sandwich::bread works", {

  # unweighted tests

  bread_lin <- bread(lin)
  bread_rob <- bread(rob)
  
  expect_equivalent(bread_lin, bread_rob) # changed from expect_equal
    
#   # weighted tests
#   
#   bread_wlm <- bread(wlin)
#   bread_wrob <- bread(wrob)
#   expect_equal(bread_wlm, bread_wrob)
#   
})


test_that("residuals_CS() works", {

  # unweighted tests

  rcs_fit <- residuals_CS(lin)
  rcs_rob <- residuals_CS(rob)

  expect_equal(rcs_fit, rcs_rob)

#   rcs_fit_fe <- residuals_CS(lin_fe)
#   rcs_rob_fe <- residuals_CS(rob_fe)
#   
#   expect_equal(rcs_fit_fe, rcs_rob_fe)
#   
#   # weighted tests
#   
#   rcs_wlm <- residuals_CS(wlin)
#   rcs_wrob <- residuals_CS(wrob)
#   
#   expect_equal(rcs_wlm, rcs_wrob)
})


test_that("coef() works", {

  # unweighted tests

  coef_fit <- coef(lin)
  coef_rob <- coef(rob)

  expect_equal(coef_fit, coef_rob)
#   
#   coef_fit_fe <- coef(lin_fe)
#   coef_rob_fe <- coef(rob_fe)
#   
#   expect_equal(coef_fit_fe[names(coef_rob_fe)], coef_rob_fe)
#   
#   # weighted tests
#   
#   coef_wlm <- coef(wlin)
#   coef_wrob <- coef(wrob)
#   
#   expect_equal(coef_wlm, coef_wrob)
#   
})


test_that("nobs() works", {

  # unweighted tests

  nobs_fit <- nobs(lin)
  nobs_rob <- nobs(rob)

  expect_equal(nobs_fit, nobs_rob)

#   # weighted tests
#   
#   nobs_wlm <- nobs(wlin)
#   nobs_wrob <- nobs(wrob)
#   
#   expect_equal(nobs_wlm, nobs_wrob)
#   
})



test_that("v_scale() works", {

  # unweighted tests

  vs_fit <- v_scale(lin)
  vs_rob <- v_scale(rob)
  
  expect_equal(vs_fit, vs_rob)
    
  # weighted tests
#   
#   vs_wlm <- v_scale(wlin)
#   vs_wrob <- v_scale(wrob)
#   
#   expect_equal(vs_wlm, vs_wrob)
#   
})


test_that("vcovCR works", {

  types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")

  # unweighted tests

  for (type in types) {

    vcov_lin <- vcovCR(lin_chick, ChickWeight$Chick, type = type)
    vcov_rob <- vcovCR(rob_chick, ChickWeight$Chick, type = type)

    expect_equal(vcov_lin, vcov_rob,
                 label = paste0("When type = ", type, ", ", "vcov_lin"))

#     
#     # weighted tests
#     
#     vcov_wlm <- vcovCR(wlin, ChickWeight$Chick, type = type)
#     vcov_wlmr <- vcovCR(wrob, ChickWeight$Chick, type = type)
#     
#     expect_equal(vcov_wlm, vcov_wlmr)
#     
#     
#     if (type %in% c("CR0","CR2")) {
#       
#       rob_type <- robust(
#         weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#         clusters = Chick, 
#         se_type = type
#       )
#       
#       expect_equal(as.matrix(vcov_lmr), vcov(rob_type), 
#                    label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr)"))
#       
#       rob_fe_type <- robust(
#         weight ~ Time:Diet, data = ChickWeight, 
#         clusters = Chick, fixed_effects = ~Chick,
#         se_type = type
#       )
#       
#       expect_equal(as.matrix(vcov_lmr_fe), vcov(rob_fe_type), 
#                    label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr_fe)"))
#       
#       wrob_type <- robust(
#         weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#         clusters = Chick, weights = wt,
#         se_type = type
#       )
#       
#       expect_equal(as.matrix(vcov_wlmr), vcov(wrob_type), 
#                    label = paste0("When type = ", type, ", ", "as.matrix(vcov_wlmr)"))
#       
#     }
  }

})


# test_that("vovCR properly pulls cluster specified for robust", {
#   
#   # unweighted tests
#   
#   uw_clust <- vcovCR(rob, ChickWeight$Chick, "CR2")
#   uw_no_clust <- vcovCR(rob, type = "CR2")
#   uw_lm <- vcovCR(lin, ChickWeight$Chick, "CR2")
#   
#   expect_equal(uw_clust, uw_no_clust)
#   expect_equal(uw_no_clust, uw_lm)
#   
#   # create an robust that draws in data differently
#   rob_fact <- robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#                            clusters = factor(ChickWeight$Chick_ordered, ordered = FALSE))
#   # perform vcovCR
#   uw_fact_cr <- vcovCR(rob_fact, type = "CR2")
#   
#   # check they are the same
#   expect_equivalent(uw_clust, uw_fact_cr)
#   
#   # put cluster data in a variable
#   # fact <- factor(ChickWeight$Chick_ordered, ordered = FALSE)
#   fact <- ChickWeight$Chick
#   
#   # pass variable to robust
#   rob_var <- robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#                           clusters = fact)
#   
#   # perform vcovCR
#   uw_fact_var <- vcovCR(rob_var, type = "CR2")
#   
#   # check they are the same
#   expect_equivalent(uw_clust, uw_fact_var)
#   
#   # weighted tests
#   
#   w_clust <- vcovCR(wrob, ChickWeight$Chick, "CR2")
#   w_no_clust <- vcovCR(wrob, type = "CR2")
#   w_lm <- vcovCR(wlin, ChickWeight$Chick, "CR2")
#   
#   expect_equal(w_clust, w_no_clust)
#   expect_equal(w_no_clust, w_lm)
#   
#   # create an robust that draws in data differently
#   rob_fact_w <- robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
#                              data = ChickWeight, 
#                              clusters = factor(ChickWeight$Chick_ordered, ordered = FALSE))
#   # perform vcovCR
#   w_fact_cr <- vcovCR(rob_fact_w, type = "CR2")
#   
#   expect_equal(w_clust, w_fact_cr)
#   
#   # pass variable to robust
#   rob_var_w <- robust(weight ~ 0 + Diet + Time:Diet, weights = wt,
#                             data = ChickWeight, clusters = fact)
#   
#   # perform vcovCR
#   w_fact_var <- vcovCR(rob_var_w, type = "CR2")
#   
#   # check they are the same
#   expect_equal(w_clust, w_fact_var)
#   
# })
# 
# 
# test_that("na.action.robust() works correctly", {
#   
#   compare_na_actions <- function(i) {
#     
#     # generate random data
#     n <- 100
#     df <- data.frame(
#       y = rnorm(n),
#       x1 = rnorm(n),
#       x2 = rnorm(n)
#     )
#     
#     # add random NA values to y and x2
#     miss_rows <- sample.int(n, size = n/10)
#     df$y[miss_rows] <- NA
#     df$x2[miss_rows] <- NA
#     
#     # fit models
#     linear <- lm(y ~ x1 + x2 + x1:x2, data = df)
#     robust <- robust(y ~ x1 + x2 + x1:x2, data = df)
#     
#     # get na.action() of models
#     na_lm <- na.action(linear)
#     na_rob <- na.action(robust)
#     
#     # compare
#     expect_equal(na_lm, na_rob)
#   }
#   
#   # compare 10 times with different random data
#   lapply(1:10, compare_na_actions)
#   
# })
# 
# 
# test_that("try_cholesky argument does not interfere with vcovCR functionality", {
#   
#   rob_chole <- robust(
#     weight ~ 0 + Diet + Time:Diet, 
#     data = ChickWeight, 
#     clusters = Chick, 
#     try_cholesky = TRUE
#   )
#   
#   expect_equal(
#     vcovCR(rob, type = "CR2"), 
#     vcovCR(rob_chole, type = "CR2")
#   )
#   
#   rob_fe_chole <- robust(
#     weight ~ Time:Diet, data = ChickWeight, 
#     clusters = Chick, 
#     fixed_effects = ~Chick,
#     try_cholesky = TRUE
#   )
#   
#   expect_equal(vcovCR(rob_fe, type = "CR2"), vcovCR(rob_fe_chole, type = "CR2"))
#   
#   wrob_chole <- robust(
#     weight ~ 0 + Diet + Time:Diet, 
#     weights = wt, 
#     data = ChickWeight, 
#     clusters = Chick,
#     try_cholesky = TRUE
#   )
#   
#   expect_equal(vcovCR(wrob, type = "CR2"), vcovCR(wrob_chole, type = "CR2"))
#   
# })


# test_that("subset argument does not interfere with vcovCR functionality", {
#   
#   lin_sub <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#                    subset = ChickWeight$rando == "Keep")
#   rob_sub <- robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
#                           clusters = Chick, subset = ChickWeight$rando == "Keep")
#   
#   expect_equal(vcovCR(lin_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2"), 
#                vcovCR(rob_sub, type = "CR2"))
#   
#   lin_fe_sub <- lm(
#     weight ~ 0 + Time:Diet + Chick, 
#     data = ChickWeight, 
#     subset = ChickWeight$rando == "Keep"
#   )
#   rob_fe_sub <- robust(
#     weight ~ Time:Diet, 
#     data = ChickWeight, 
#     clusters = Chick, 
#     fixed_effects = ~Chick,
#     subset = ChickWeight$rando == "Keep"
#   )
#   
#   sub_coef <- names(coef(rob_fe_sub))
#   expect_equivalent(
#     vcovCR(lin_fe_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2")[sub_coef, sub_coef],
#     as.matrix(vcovCR(rob_fe_sub, type = "CR2"))
#   )
#   
#   wlin_sub <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, 
#                     data = ChickWeight, subset = ChickWeight$rando == "Keep")
#   wrob_sub <- robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
#                            data = ChickWeight, clusters = Chick, 
#                            subset = ChickWeight$rando == "Keep")
#   
#   expect_equal(vcovCR(wlin_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2"),
#                vcovCR(wrob_sub, type = "CR2"))
#   
# })


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

# Create identical lm and robust models
belts_fit <- lm(DriversKilled ~ kms + PetrolPrice + law + year, data = belts)
belts_rob <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts)


test_that("Wald_test() works with robust", {
  
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


test_that("conf_int() works with robust", {
  
  conf_FIT <- conf_int(belts_fit, vcov = "CR2", cluster = belts$month)
  conf_ROB <- conf_int(belts_rob, vcov = "CR2", cluster = belts$month)
  
  expect_equal(conf_FIT, conf_ROB)
  
})


test_that("coef_test() works with robust", {
  
  coef_FIT <- coef_test(belts_fit, vcov = "CR2", cluster = belts$month)
  coef_ROB <- coef_test(belts_rob, vcov = "CR2", cluster = belts$month)
  
  expect_equal(coef_FIT, coef_ROB)
  
})


# =============== Tests Based on test_lm.R ===============


test_that("Order doesn't matter.",{
  
  check_sort_order(belts_rob, belts, "month", tol = 1e-5)
  
})


test_that("clubSandwich works with dropped observations", {
  belts_miss <- belts
  miss_indicator <- sample.int(nrow(belts), size = round(nrow(belts) / 10))
  belts_miss$law[miss_indicator] <- NA
  belts_miss$kms[miss_indicator] <- NA
  
  rob_dropped <- lm_robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts_miss, cluster = month)
  belts_complete <- subset(belts_miss, !is.na(law))
  rob_complete <- robust(DriversKilled ~ kms + PetrolPrice + law + year, data = belts_complete, cluster = month)
  
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
