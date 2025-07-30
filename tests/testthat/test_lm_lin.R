context("lm_lin robust objects")

skip_if_not_installed("estimatr")

suppressMessages(library(estimatr, quietly=TRUE))

set.seed(20250630)
data("CO2")
N <- nrow(CO2)
CO2$wt <- 1 + rpois(N, 3)
CO2$Plant_[!is.na(CO2$uptake), ] # remove rows with NA in uptakeordered <- CO2$Plant
CO2$Plant <- factor(CO2$Plant, ordered = FALSE)
CO2$Plant_int <- as.integer(CO2$Plant)
CO2$rando <- "Drop"
CO2$rando[sample(1:N, size = round(0.8 * N))] <- "Keep"
table(CO2$rando)

# make a copy of CO2 for weighted clustered models
wcCO2 <- CO2

CO2$conc_c <- CO2$conc - mean(CO2$conc)
CO2$`(factor(Type)Mississippi)_c` <- ifelse(CO2$Type == "Mississippi", 0.5, -0.5)

lin <- lm_lin(
  uptake ~ Treatment,
  covariates = ~ conc,
  data = CO2
)

lin_exp <- lm_lin(
  uptake ~ Treatment,
  covariates = ~ conc + factor(Type),
  data = CO2
)

rob <- lm_robust(
  uptake ~ Treatment + 
    conc_c + Treatment:conc_c,
  data = CO2
)

rob_exp <- lm_robust(
  uptake ~ Treatment + conc_c + `(factor(Type)Mississippi)_c` + Treatment:conc_c + Treatment:`(factor(Type)Mississippi)_c`,
  data = CO2
)

# weighted and clustered

wclin <- lm_lin(
  uptake ~ Treatment,
  covariates = ~ conc,
  data = wcCO2,
  weights = wt,
  clusters = Plant
)

wcCO2$conc_c <- wcCO2$conc - wclin$scaled_center[["conc"]]

wcrob <- lm_robust(
  uptake ~ Treatment +
    conc_c + Treatment:conc_c,
  data = wcCO2,
  weights = wt,
  clusters = Plant
)

test_that("model.frame() works", {

  # basic models
  mf_lin <- model.frame(lin)
  mf_rob <- model.frame(rob)
  expect_equivalent(mf_lin, mf_rob)
  
  # covariates with an expression 
  mf_exp <- model.frame(lin_exp)
  mf_rexp <- model.frame(rob_exp)
  expect_equivalent(mf_exp, mf_rexp)
  
  # weighted models
  mf_wclin <- model.frame(wclin)
  mf_wcrob <- model.frame(wcrob)
  expect_equivalent(mf_wclin, mf_wcrob[names(mf_wclin)])

})

test_that("model_matrix() works", {
  
  # basic models
  mm_lin <- model_matrix(lin)
  mm_rob <- model_matrix(rob)
  expect_equivalent(mm_lin, mm_rob)
  
  # covariates with an expression 
  mm_exp <- model_matrix(lin_exp)
  mm_rexp <- model_matrix(rob_exp)
  expect_equivalent(mm_exp, mm_rexp)
  
  # weighted models
  mm_wclin <- model_matrix(wclin)
  mm_wcrob <- model_matrix(wcrob)
  expect_equivalent(mm_wclin, mm_wcrob)
  
})

test_that("targetVariance() works", {

  # unweighted tests

  tV_lin <- targetVariance(lin, CO2$Plant)
  tV_rob <- targetVariance(rob, CO2$Plant)

  expect_equal(tV_lin, tV_rob)
  
  # weighted tests
  
  tV_wclin <- targetVariance(wclin, wcCO2$Plant)
  tV_wcrob <- targetVariance(wcrob, wcCO2$Plant)
  
  expect_equal(tV_wclin, tV_wcrob)
  
})


test_that("weightMatrix() works", {

  # unweighted tests

  wM_fit <- weightMatrix(lin, CO2$Plant)
  wM_rob <- weightMatrix(rob, CO2$Plant)

  expect_equal(wM_fit, wM_rob)

  # weighted tests

  wM_wcfit <- weightMatrix(wclin, wcCO2$Plant)
  wM_wcrob <- weightMatrix(wcrob, wcCO2$Plant)
  
  expect_equal(wM_wcfit, wM_wcrob)
  
})

test_that("sandwich::bread works", {

  # unweighted tests

  bread_lin <- bread(lin)
  bread_rob <- bread(rob)
  
  expect_equivalent(bread_lin, bread_rob) # changed from expect_equal
    
  # weighted tests
  
  bread_wclin <- bread(wclin)
  bread_wcrob <- bread(wcrob)
  
  expect_equivalent(bread_wclin, bread_wcrob) # changed from expect_equal

})


test_that("residuals_CS() works", {

  # unweighted tests

  rcs_lin <- residuals_CS(lin)
  rcs_rob <- residuals_CS(rob)

  expect_equal(rcs_lin, rcs_rob)
 
  # weighted tests
  
  rcs_wclin <- residuals_CS(wclin)
  rcs_wcrob <- residuals_CS(wcrob)
  
  expect_equal(rcs_wclin, rcs_wcrob)

})


test_that("coef() works", {

  # unweighted tests

  coef_lin <- coef(lin)
  coef_rob <- coef(rob)

  expect_equal(coef_lin, coef_rob)
    
  # weighted tests

  coef_wclin <- coef(wclin)
  coef_wcrob <- coef(wcrob)
  
  expect_equal(coef_wclin, coef_wcrob)
    
})


test_that("nobs() works", {

  # unweighted tests

  nobs_lin <- nobs(lin)
  nobs_rob <- nobs(rob)

  expect_equal(nobs_lin, nobs_rob)

  # weighted tests

  nobs_wclin <- nobs(wclin)
  nobs_wcrob <- nobs(wcrob)
  
  expect_equal(nobs_wclin, nobs_wcrob)
  
})



test_that("v_scale() works", {

  # unweighted tests

  vs_lin <- v_scale(lin)
  vs_rob <- v_scale(rob)
  
  expect_equal(vs_lin, vs_rob)
    
  # weighted tests
  
  vs_wclin <- v_scale(wclin)
  vs_wcrob <- v_scale(wcrob)
  
  expect_equal(vs_wclin, vs_wcrob)
  
})


test_that("vcovCR works", {

  types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")

  for (type in types) {

    # unweighted tests
    
    vcov_lin <- vcovCR(lin, CO2$Plant, type = type)
    vcov_rob <- vcovCR(rob, CO2$Plant, type = type)

    # changed from equal b/c of attributes
    expect_equivalent(vcov_lin, vcov_rob,
                 label = paste0("When type = ", type, ", ", "vcov_lin"))

    # weighted tests

    vcov_wclin <- vcovCR(wclin, wcCO2$Plant, type = type)
    vcov_wcrob <- vcovCR(wcrob, wcCO2$Plant, type = type)

    # changed from equal b/c of attributes
    expect_equivalent(vcov_wclin, vcov_wcrob,
                 label = paste0("When type = ", type, ", ", "vcov_wclin"))

  }

})


test_that("vovCR properly pulls cluster specified for lm_lin objects", {

  # weighted tests

  w_clust <- vcovCR(wclin, wcCO2$Plant, "CR2")
  w_no_clust <- vcovCR(wclin, type = "CR2")
  w_lin <- vcovCR(wclin, wcCO2$Plant, "CR2")

  expect_equal(w_clust, w_no_clust)
  expect_equal(w_no_clust, w_lin)

  # create an lm_robust object that draws in data differently
  lin_fact <- lm_lin(uptake ~ Treatment,
                     covariates = ~ conc,
                     data = CO2,
                     weights = wt,
                     clusters = factor(CO2$Plant, ordered = FALSE))
  
  # perform vcovCR
  w_fact_cr <- vcovCR(lin_fact, type = "CR2")

  # check they are the same
  expect_equivalent(w_clust, w_fact_cr)

  # put cluster data in a variable
  # fact <- factor(ChickWeight$Chick_ordered, ordered = FALSE)
  fact <- CO2$Plant

  # pass variable to robust
  lin_var <- lm_lin(uptake ~ Treatment,
                    covariates = ~ conc,
                    data = CO2,
                    weights = wt,
                    clusters = fact)

  # perform vcovCR
  w_fact_var <- vcovCR(lin_var, type = "CR2")

  # check they are the same
  expect_equivalent(w_clust, w_fact_var)

})


test_that("vcovCR works with se_type inherited from lm_lin().", {
  
  types <- c("CR0", "CR2", "stata")
  
  data("AchievementAwardsRCT")
  AchievementAwardsRCT$sub <- with(
    AchievementAwardsRCT, 
    unsplit(tapply(Bagrut_status, school_id, \(x) seq(0, 1, length.out = length(x))), school_id)
  )
  AchievementAwardsRCT <- subset(AchievementAwardsRCT, school_id != 1 & sub < 0.2)
  AchievementAwardsRCT$father_ed <- ifelse(AchievementAwardsRCT$father_ed <= 5, "GS", ifelse(AchievementAwardsRCT$father_ed <=12, "HS","Coll"))
  
  for (type in types) {
    
    # basic
    rob_fit <- lm_lin(
      achv_math ~ treated,
      covariates = ~ school_type + sex + father_ed,
      data = AchievementAwardsRCT, 
      cluster = school_id,
      se_type = type
    )
    
    expect_equal(as.matrix(vcovCR(rob_fit)), vcov(rob_fit))
    
    # with expression in the covariates
    rob_exp <- lm_lin(
      achv_math ~ treated,
      covariates = ~ school_type + sex + factor(father_ed),
      data = AchievementAwardsRCT, 
      cluster = school_id,
      se_type = type
    )
    
    expect_equivalent(as.matrix(vcovCR(rob_exp)), vcov(rob_exp))
    
    # weighted
    wt_fit <- lm_lin(
      achv_math ~ treated,
      covariates = ~ school_type + sex + father_ed,
      data = AchievementAwardsRCT, 
      weights = siblings + 1,
      cluster = school_id,
      se_type = type
    )
    
    expect_equal(as.matrix(vcovCR(wt_fit)), vcov(wt_fit))
    
  }
  
})

test_that("na.action.robust() works correctly", {

  compare_na_actions <- function(i) {

    # generate random data
    n <- 100
    df <- data.frame(
      y = rnorm(n),
      x1 = sample(letters[1:2], size = n, replace = TRUE),
      x2 = rnorm(n)
    )

    # add random NA values to y and x2
    miss_rows <- sample.int(n, size = n/10)
    df$y[miss_rows] <- NA
    df$x2[miss_rows] <- NA

    # fit models
    lin_fit <- lm_lin(y ~ x1,
                  covariates = ~ x2,
                  data = df)
    rob_fit <- lm_robust(y ~ x1 + x2 + x1:x2,
                     data = df)

    # get na.action() of models
    na_lin <- na.action(lin_fit)
    na_rob <- names(na.action(rob_fit))

    # compare
    expect_equal(na_lin, na_rob)
  }

  # compare 10 times with different random data
  lapply(1:10, compare_na_actions)

})


test_that("try_cholesky argument does not interfere with vcovCR functionality", {

  lin_chole <- lm_lin(
    uptake ~ Treatment,
    covariates = ~ conc,
    data = CO2,
    try_cholesky = TRUE
  )

  expect_equal(
    vcovCR(lin, cluster = CO2$Plant, type = "CR2"),
    vcovCR(lin_chole, cluster = CO2$Plant, type = "CR2")
  )

  wclin_chole <- wclin <- lm_lin(
    uptake ~ Treatment,
    covariates = ~ conc,
    data = wcCO2,
    weights = wt,
    clusters = Plant,
    try_cholesky = TRUE
  )

  expect_equal(vcovCR(wclin, type = "CR2"), vcovCR(wclin_chole, type = "CR2"))

})


test_that("subset argument does not interfere with vcovCR functionality", {

  lin_sub <- lm_lin(uptake ~ Treatment,
                    covariates = ~ conc,
                    data = CO2,
                    subset = rando == "Keep")
  
  # CO2$conc_c_no_subset <- CO2$conc_c
  CO2$conc_c <- CO2$conc - lin_sub$scaled_center
  
  rob_sub <- lm_robust(uptake ~ Treatment + conc_c + Treatment:conc_c,
                       data = CO2,
                       subset = rando == "Keep")

  # changed from equal b/c of attributes
  expect_equivalent(vcovCR(lin_sub, CO2$Plant[CO2$rando == "Keep"], type = "CR2"),
               vcovCR(rob_sub, CO2$Plant[CO2$rando == "Keep"], type = "CR2"))

  wclin_sub <- lm_lin(uptake ~ Treatment,
                      covariates = ~ conc,
                      data = wcCO2,
                      weights = wt,
                      clusters = Plant,
                      subset = rando == "Keep")
  
  # wcCO2$conc_c_no_subset <- wcCO2$conc_c
  wcCO2$conc_c <- wcCO2$conc - wclin_sub$scaled_center
  
  wcrob_sub <- lm_robust(uptake ~ Treatment + conc_c + Treatment:conc_c,
                      data = wcCO2,
                      weights = wt,
                      clusters = Plant,
                      subset = rando == "Keep")

  # changed from equal b/c of attributes
  expect_equivalent(vcovCR(wclin_sub, type = "CR2"),
               vcovCR(wcrob_sub, type = "CR2"))

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
