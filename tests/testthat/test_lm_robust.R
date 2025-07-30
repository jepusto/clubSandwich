context("lm_robust objects")

skip_if_not_installed("estimatr")

suppressMessages(library(estimatr, quietly=TRUE))


set.seed(20250630)
data("ChickWeight", package = "datasets")
N <- nrow(ChickWeight)
ChickWeight$wt <- 1 + rpois(N, 3)
ChickWeight$Chick_ordered <- ChickWeight$Chick # James' suggestion 4/16
ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)
ChickWeight$Chick_int <- as.integer(ChickWeight$Chick)
ChickWeight$rando <- "Drop"
ChickWeight$rando[sample(1:N, size = round(0.8 * N))] <- "Keep"
table(ChickWeight$rando)

lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(
  weight ~ 0 + Diet + Time:Diet, 
  data = ChickWeight, 
  clusters = Chick
)

lm_fit_fe <- lm(weight ~ 0 + Time:Diet + Chick, data = ChickWeight)
lm_rob_fe <- lm_robust(
  weight ~ Time:Diet, 
  data = ChickWeight, 
  clusters = Chick, 
  fixed_effects = ~ Chick
)

wlm_fit <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, data = ChickWeight)
wlm_rob <- lm_robust(
  weight ~ 0 + Diet + Time:Diet, 
  weights = wt, 
  data = ChickWeight, 
  clusters = Chick
)

# Two-way fixed effects models
data("MortalityRates")
MortalityRates <- subset(MortalityRates, cause == "Motor Vehicle")
MortalityRates$state <- factor(MortalityRates$state)
MortalityRates$year <- factor(MortalityRates$year)

lm_fit_fe2 <- lm(
  mrate ~ 0 + year + state + legal + beertaxa + beerpercap + winepercap + spiritpercap,
  data = MortalityRates
)
lm_rob_fe2 <- lm_robust(
  mrate ~ legal + beertaxa + beerpercap + winepercap + spiritpercap,
  fixed_effects = ~ year + state,
  data = MortalityRates,
  cluster = state
)

test_that("model.frame() works", {
  
  # unweighted models
  mf_fit <- model.frame(lm_fit)
  mf_rob <- model.frame(lm_rob)
  mf_rob$`(clusters)` <- NULL
  expect_equal(mf_fit, mf_rob)

  # fixed effects models
  mf_fit_fe <- model.frame(lm_fit_fe)
  mf_rob_fe <- model.frame(lm_rob_fe)
  mf_rob_fe$`(clusters)` <- NULL
  expect_equivalent(mf_fit_fe, mf_rob_fe)

  # two-way fixed effects models
  mf_fit_fe2 <- model.frame(lm_fit_fe2)
  mf_rob_fe2 <- model.frame(lm_rob_fe2)
  mf_rob_fe2$`(clusters)` <- NULL
  expect_equivalent(mf_fit_fe2[names(mf_rob_fe2)], mf_rob_fe2)
  
  # weighted models
  mf_wlm <- model.frame(wlm_fit)
  mf_wrob <- model.frame(wlm_rob)
  mf_wrob$`(clusters)` <- NULL

  expect_equal(mf_wlm, mf_wrob)

  # set up tests with missing values
  set.seed(20250629)
  
  dat_miss <- ChickWeight
  i <- 1:nrow(ChickWeight)
  miss1 <- 3:8
  miss2 <- sample(i, 9L)
  dat_miss$Diet_miss <- dat_miss$Diet
  dat_miss$Diet_miss[miss1] <- NA
  dat_miss$wt_miss <- ifelse(i %in% miss1, NA, dat_miss$wt)
  dat_miss$Chick_miss1 <- dat_miss$Chick
  dat_miss$Chick_miss1[miss1] <- NA
  dat_miss$Chick_miss2 <- dat_miss$Chick
  dat_miss$Chick_miss2[miss2] <- NA
  
  dat_complete1 <- droplevels(dat_miss[setdiff(i, miss1),])
  dat_complete2 <- droplevels(dat_miss[setdiff(i, miss2),])
  dat_complete <- droplevels(dat_miss[setdiff(i, c(miss1, miss2)),])
  
  compare_model_frames <- function(data1, data2, ...) {
    cl <- match.call()
    
    cl1 <- cl
    cl1$data2 <- NULL
    names(cl1)[2] <- "data"
    cl1[[1]] <- quote(estimatr::lm_robust)
    m1 <- suppressWarnings(eval(cl1, parent.frame()))
    mf1 <- model.frame(m1)
    
    cl2 <- cl
    cl2$data1 <- NULL
    names(cl2)[2] <- "data"
    cl2[[1]] <- quote(estimatr::lm_robust)
    m2 <- suppressWarnings(eval(cl2, parent.frame()))
    mf2 <- model.frame(m2)
    
    expect_equivalent(mf1, mf2)
  }
  
  # X missing
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet_miss, 
    fixed_effects = ~ Chick + Diet,
    clusters = Chick_int,
    se_type = "CR0"
  )
  
  # clusters missing
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet, 
    fixed_effects = ~ Chick,
    clusters = Chick_miss1,
    se_type = "CR0"
  )
  
  # weights missing  
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet, 
    weights = wt_miss,
    fixed_effects = ~ Chick,
    se_type = "HC0"
  )
  
  # FE missing
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet, 
    fixed_effects = ~ Chick + Diet_miss,
    clusters = Chick_int,
    se_type = "CR0"
  )
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete2, 
    formula = weight ~ Time:Diet, 
    fixed_effects = ~ Chick_miss2 + Diet,
    se_type = "HC0"
  )
  
  # X and FE missing
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet_miss, 
    fixed_effects = ~ Chick_miss1,
    clusters = Chick_int,
    se_type = "CR0"
  )
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete, 
    formula = weight ~ Time:Diet_miss, 
    fixed_effects = ~ Chick_miss2,
    clusters = Chick_int,
    se_type = "CR0"
  )
  
  # clusters and FE missing
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet, 
    fixed_effects = ~ Chick_miss1,
    clusters = Chick_miss1,
    se_type = "CR0"
  )
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete, 
    formula = weight ~ Time:Diet, 
    fixed_effects = ~ Chick_miss1,
    clusters = Chick_miss2,
    se_type = "CR0"
  )
  
  # weights and FE missing  
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete1, 
    formula = weight ~ Time:Diet, 
    weights = wt_miss,
    fixed_effects = ~ Chick_miss1,
    clusters = Chick_int,
    se_type = "CR0"
  )
  compare_model_frames(
    data1 = dat_miss, data2 = dat_complete, 
    formula = weight ~ Time:Diet, 
    weights = wt_miss,
    fixed_effects = ~ Chick_miss2,
    se_type = "HC0"
  )
})

test_that("model_matrix() works", {
  
  # basic models
  
  mm_fit <- model_matrix(lm_fit) 
  mm_rob <- model_matrix(lm_rob)
  expect_equivalent(mm_fit, mm_rob)

  
  # fixed effects models
  
  mm_fit_fe <- model_matrix(lm_fit_fe)
  mm_rob_fe <- model_matrix(lm_rob_fe)
  amm_rob_fe <- augmented_model_matrix(lm_rob_fe)

  # Check that fixed effects are the same
  expect_equivalent(mm_fit_fe[,colnames(amm_rob_fe)], amm_rob_fe)
  # Core predictor matrices are different
  expect_false(identical(mm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  # But one can be computed by residualizing
  expect_equivalent(
    mm_rob_fe,
    residuals(lm.fit(amm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  )
  # model matrix is centered by chick so means are zero
  expect_lt(max(abs(apply(mm_rob_fe, 2, \(x) tapply(x, ChickWeight$Chick, mean)))), 1e-12)
  
  
  # two-way fixed effects models
  
  mm_fit_fe2 <- model_matrix(lm_fit_fe2)
  mm_rob_fe2 <- model_matrix(lm_rob_fe2)
  amm_rob_fe2 <- augmented_model_matrix(lm_rob_fe2)
  
  # Check that fixed effects are the same
  expect_equivalent(mm_fit_fe2[,colnames(amm_rob_fe2)], amm_rob_fe2)
  
  # Core predictor matrices are different
  expect_false(identical(mm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  # But one can be computed by residualizing
  expect_equivalent(
    mm_rob_fe2,
    residuals(lm.fit(amm_rob_fe2, mm_fit_fe2[,colnames(mm_rob_fe2)]))
  )
  # model matrix is centered by state and year
  MortalityRates_full <- subset(MortalityRates, !is.na(beertaxa))
  expect_lt(max(abs(apply(mm_rob_fe2, 2, \(x) tapply(x, MortalityRates_full$state, mean)))), 1e-12)
  expect_lt(max(abs(apply(mm_rob_fe2, 2, \(x) tapply(x, MortalityRates_full$year, mean)))), 1e-12)

    
  # weighted models
  
  mm_wlm <- model_matrix(wlm_fit)
  mm_wrob <- model_matrix(wlm_rob)
  expect_equivalent(mm_wlm, mm_wrob)
  
})

test_that("model_matrix() works without fixest", {
  
  local_mocked_bindings(
    requireNamespace = function(...) FALSE
  )
  
  expect_false(requireNamespace("fixest", quietly = TRUE))
  
  # basic models
  
  mm_fit <- model_matrix(lm_fit) 
  mm_rob <- model_matrix(lm_rob)
  amm_rob_fe <- augmented_model_matrix(lm_rob_fe)
  expect_equivalent(mm_fit, mm_rob)
  
  
  # fixed effects models
  
  mm_fit_fe <- model_matrix(lm_fit_fe)
  mm_rob_fe <- model_matrix(lm_rob_fe)
  
  # Check that fixed effects are the same
  expect_equivalent(mm_fit_fe[,colnames(amm_rob_fe)], amm_rob_fe)
  # Core predictor matrices are different
  expect_false(identical(mm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  # But one can be computed by residualizing
  expect_equivalent(
    mm_rob_fe,
    residuals(lm.fit(amm_rob_fe, mm_fit_fe[,colnames(mm_rob_fe)]))
  )
  # model matrix is centered by chick so means are zero
  expect_lt(max(abs(apply(mm_rob_fe, 2, \(x) tapply(x, ChickWeight$Chick, mean)))), 1e-12)
  
  expect_message(
    vcovCR(lm_rob_fe,type = "CR0"),
    "For improved performance in models with fixed effects, install the package \\{fixest\\}\\."
  )
  
})

test_that("targetVariance() works", {
  
  # unweighted tests
  
  tV_fit <- targetVariance(lm_fit, ChickWeight$Chick)
  tV_rob <- targetVariance(lm_rob, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob)
  
  tV_rob_fe <- targetVariance(lm_rob_fe, ChickWeight$Chick)
  
  expect_equal(tV_fit, tV_rob_fe)
  
  # weighted tests
  
  tV_wlm <- targetVariance(wlm_fit, ChickWeight$Chick)
  tV_wrob <- targetVariance(wlm_rob, ChickWeight$Chick)
  
  expect_equal(tV_wlm, tV_wrob)
  
})


test_that("weightMatrix() works", {
  
  # unweighted tests
  
  wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick)
  wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob)
  
  wM_rob_fe <- weightMatrix(lm_rob_fe, ChickWeight$Chick)
  
  expect_equal(wM_fit, wM_rob_fe)
  
  # weighted tests
  
  wM_wlm <- weightMatrix(wlm_fit, ChickWeight$Chick)
  wM_wrob <- weightMatrix(wlm_rob, ChickWeight$Chick)
  
  expect_equal(wM_wlm, wM_wrob)
  
})

test_that("sandwich::bread works", {
  
  # basic models
  bread_lm <- bread(lm_fit)
  bread_rob <- bread(lm_rob)
  expect_equal(bread_lm, bread_rob)
  
  # fixed effects models
  bread_lm_fe <- bread(lm_fit_fe)
  bread_rob_fe <- bread(lm_rob_fe)
  focal_coefs <- names(coef(lm_rob_fe))
  expect_equal(bread_lm_fe[focal_coefs,focal_coefs], as.matrix(bread_rob_fe))
  
  # two-way fixed effects models
  bread_lm_fe2 <- bread(lm_fit_fe2)
  bread_rob_fe2 <- bread(lm_rob_fe2)
  focal_coefs <- names(coef(lm_rob_fe2))
  expect_equal(bread_lm_fe2[focal_coefs,focal_coefs], as.matrix(bread_rob_fe2))
  
  # weighted models
  bread_wlm <- bread(wlm_fit)
  bread_wrob <- bread(wlm_rob)
  expect_equal(bread_wlm, bread_wrob)
  
})


test_that("residuals_CS() works", {
  
  # basic models
  rcs_fit <- residuals_CS(lm_fit)
  rcs_rob <- residuals_CS(lm_rob)
  expect_equal(rcs_fit, rcs_rob)
  
  # fixed effects models
  rcs_fit_fe <- residuals_CS(lm_fit_fe)
  rcs_rob_fe <- residuals_CS(lm_rob_fe)
  expect_equal(rcs_fit_fe, rcs_rob_fe)
  
  # two-way fixed effects models
  rcs_fit_fe2 <- residuals_CS(lm_fit_fe2)
  rcs_rob_fe2 <- residuals_CS(lm_rob_fe2)
  expect_equal(rcs_fit_fe2, rcs_rob_fe2)

  # weighted models
  rcs_wlm <- residuals_CS(wlm_fit)
  rcs_wrob <- residuals_CS(wlm_rob)
  expect_equal(rcs_wlm, rcs_wrob)
  
})


test_that("coef() works", {
  
  # basic models
  coef_fit <- coef(lm_fit)
  coef_rob <- coef(lm_rob)
  expect_equal(coef_fit, coef_rob)
  
  # fixed effects models
  coef_fit_fe <- coef(lm_fit_fe)
  coef_rob_fe <- coef(lm_rob_fe)
  expect_equal(coef_fit_fe[names(coef_rob_fe)], coef_rob_fe)
  
  # two-way fixed effects models
  coef_fit_fe2 <- coef(lm_fit_fe2)
  coef_rob_fe2 <- coef(lm_rob_fe2)
  expect_equal(coef_fit_fe2[names(coef_rob_fe2)], coef_rob_fe2)
  
  # weighted models
  coef_wlm <- coef(wlm_fit)
  coef_wrob <- coef(wlm_rob)
  expect_equal(coef_wlm, coef_wrob)
  
})


test_that("nobs() works", {
  
  # unweighted tests
  
  nobs_fit <- nobs(lm_fit)
  nobs_rob <- nobs(lm_rob)
  
  expect_equal(nobs_fit, nobs_rob)
  
  nobs_rob_fe <- nobs(lm_rob_fe)
  
  expect_equal(nobs_fit, nobs_rob_fe)
  
  # weighted tests
  
  nobs_wlm <- nobs(wlm_fit)
  nobs_wrob <- nobs(wlm_rob)
  
  expect_equal(nobs_wlm, nobs_wrob)
  
})



test_that("v_scale() works", {
  
  # unweighted tests
  
  vs_fit <- v_scale(lm_fit)
  vs_rob <- v_scale(lm_rob)
  vs_rob_fe <- v_scale(lm_rob_fe)
  
  expect_equal(vs_fit, vs_rob)
  expect_equal(vs_fit, vs_rob_fe)
  
  # weighted tests
  
  vs_wlm <- v_scale(wlm_fit)
  vs_wrob <- v_scale(wlm_rob)
  
  expect_equal(vs_wlm, vs_wrob)
  
})


test_that("vcovCR works", {
  
  types <- c("CR0", "CR1", "CR1p", "CR1S", "CR2", "CR3")
  
  
  focal_coefs <- names(coef(lm_rob_fe))
  
  for (type in types) {
    
    # basic models
    vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, type = type)
    vcov_lmr <- vcovCR(lm_rob, ChickWeight$Chick, type = type)
    expect_equal(vcov_lm, vcov_lmr, 
                 label = paste0("When type = ", type, ", ", "vcov_lm"))
    
    
    # fixed effects models
    
    if (type %in% c("CR0", "CR1", "CR2")) {
      
      vcov_lm_fe <- vcovCR(lm_fit_fe, ChickWeight$Chick, type = type)
      vcov_lmr_fe <- vcovCR(lm_rob_fe, ChickWeight$Chick, type = type)
      
      expect_equal(vcov_lm_fe[focal_coefs,focal_coefs], as.matrix(vcov_lmr_fe), 
                   label = paste0("When type = ", type, ", ", "vcov_lm_fe[focal_coefs,focal_coefs]"))
      
    }
    
    # weighted models
    vcov_wlm <- vcovCR(wlm_fit, ChickWeight$Chick, type = type)
    vcov_wlmr <- vcovCR(wlm_rob, ChickWeight$Chick, type = type)
    expect_equal(vcov_wlm, vcov_wlmr)
    
    
    if (type %in% c("CR0","CR2")) {
      
      # basic models
      lm_rob_type <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, 
        se_type = type
      )
      
      expect_equal(as.matrix(vcov_lmr), vcov(lm_rob_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr)"))
      
      # fixed effects models
      lm_rob_fe_type <- lm_robust(
        weight ~ Time:Diet, data = ChickWeight, 
        clusters = Chick, fixed_effects = ~Chick,
        se_type = type
      )
      
      expect_equal(as.matrix(vcov_lmr_fe), vcov(lm_rob_fe_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_lmr_fe)"))
      

      # weighted models
      
      wlm_rob_type <- lm_robust(
        weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
        clusters = Chick, weights = wt,
        se_type = type
      )
      
      expect_equal(as.matrix(vcov_wlmr), vcov(wlm_rob_type), 
                   label = paste0("When type = ", type, ", ", "as.matrix(vcov_wlmr)"))
      
    }
  }
  
})


test_that("vcovCR works with se_type inherited from lm_robust().", {
  
  types <- c("CR0", "CR2", "stata")
  
  data("OrchardSprays", package = "datasets")
  OrchardSprays$wt <- 1 + rpois(nrow(OrchardSprays), lambda = 3)

  for (type in types) {
    
    # unweighted
    rob_fit <- lm_robust(
      decrease ~ 0 + factor(rowpos) + treatment, 
      data = OrchardSprays, 
      cluster = colpos,
      se_type = type
    )
    
    expect_equal(as.matrix(vcovCR(rob_fit)), vcov(rob_fit))

    # weighted
    wt_fit <- lm_robust(
      decrease ~ 0 + factor(rowpos) + treatment, 
      weights = wt,
      data = OrchardSprays, 
      cluster = colpos,
      se_type = type
    )
    
    expect_equal(as.matrix(vcovCR(wt_fit)), vcov(wt_fit))
    
    
    if (type != "stata") {
      
      # fixed effects
      fe_fit <- lm_robust(
        decrease ~ treatment, 
        fixed_effects = ~ rowpos,
        data = OrchardSprays, 
        cluster = colpos,
        se_type = type
      )
      expect_equal(as.matrix(vcovCR(fe_fit)), vcov(fe_fit))
      
      # two-way fixed effects
      fe2_fit <- lm_robust(
        mrate ~ legal + beertaxa + beerpercap + winepercap + spiritpercap,
        fixed_effects = ~ year + state,
        data = MortalityRates,
        cluster = state,
        se_type = type
      )
      expect_equal(as.matrix(vcovCR(fe2_fit)), vcov(fe2_fit))
      
    }
    
    
  }
  
})


test_that("vovCR properly pulls cluster specified for lm_robust", {
  
  # basic models
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


test_that("try_cholesky argument does not interfere with vcovCR functionality", {
  
  lm_rob_chole <- lm_robust(
    weight ~ 0 + Diet + Time:Diet, 
    data = ChickWeight, 
    clusters = Chick, 
    try_cholesky = TRUE
  )
  
  expect_equal(
    vcovCR(lm_rob, type = "CR2"), 
    vcovCR(lm_rob_chole, type = "CR2")
  )
  
  lm_rob_fe_chole <- lm_robust(
    weight ~ Time:Diet, data = ChickWeight, 
    clusters = Chick, 
    fixed_effects = ~Chick,
    try_cholesky = TRUE
  )
  
  expect_equal(vcovCR(lm_rob_fe, type = "CR2"), vcovCR(lm_rob_fe_chole, type = "CR2"))
  
  wlm_rob_chole <- lm_robust(
    weight ~ 0 + Diet + Time:Diet, 
    weights = wt, 
    data = ChickWeight, 
    clusters = Chick,
    try_cholesky = TRUE
  )
  
  expect_equal(vcovCR(wlm_rob, type = "CR2"), vcovCR(wlm_rob_chole, type = "CR2"))
  
})


test_that("subset argument does not interfere with vcovCR functionality", {
  
  # basic models
  lm_fit_sub <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                            subset = ChickWeight$rando == "Keep")
  lm_rob_sub <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight, 
                            clusters = Chick, subset = ChickWeight$rando == "Keep")
  
  expect_equal(vcovCR(lm_fit_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2"), 
               vcovCR(lm_rob_sub, type = "CR2"))
  
  
  # fixed effects models
  
  lm_fit_fe_sub <- lm(
    weight ~ 0 + Time:Diet + Chick, 
    data = ChickWeight, 
    subset = rando == "Keep"
  )
  lm_rob_fe_sub <- lm_robust(
    weight ~ Time:Diet, 
    data = ChickWeight, 
    clusters = Chick, 
    fixed_effects = ~Chick,
    subset = rando == "Keep"
  )
  
  sub_coef <- names(coef(lm_rob_fe_sub))
  expect_equivalent(
    vcovCR(lm_fit_fe_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2")[sub_coef, sub_coef],
    as.matrix(vcovCR(lm_rob_fe_sub, type = "CR2"))
  )
  
  # weighted models
  
  wlm_fit_sub <- lm(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                            data = ChickWeight, subset = rando == "Keep")
  wlm_rob_sub <- lm_robust(weight ~ 0 + Diet + Time:Diet, weights = wt, 
                            data = ChickWeight, clusters = Chick, 
                            subset = rando == "Keep")
  
  expect_equal(vcovCR(wlm_fit_sub, ChickWeight$Chick[ChickWeight$rando == "Keep"], type = "CR2"),
               vcovCR(wlm_rob_sub, type = "CR2"))
  
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


test_that("Wald_test() works with lm_robust", {
  
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


test_that("conf_int() works with lm_robust", {
  
  conf_FIT <- conf_int(belts_fit, vcov = "CR2", cluster = belts$month)
  conf_ROB <- conf_int(belts_rob, vcov = "CR2", cluster = belts$month)
  
  expect_equal(conf_FIT, conf_ROB)
  
})


test_that("coef_test() works with lm_robust", {
  
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
