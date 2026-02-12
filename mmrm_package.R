library(mmrm)

# ===============================================
# 1: Data used in the official package website
# ===============================================

# fev_data is a clinical trial dataset
data(fev_data)

str(fev_data)
head(fev_data, 20)
# 800 total rows with 10 columns

# repeated measures means multiple rows per subject
table(fev_data$USUBJID)[1:10]  # 4 rows per subject
length(unique(fev_data$USUBJID))  # 200 subjects total
table(fev_data$AVISIT)  # 200 observations/rows per visit (4 Visits)

# Check if there are some NA's in FEV1 (forced expired volume in one second) is a measure of how quickly the lungs can be emptied
sum(is.na(fev_data$FEV1))

# mmrm handles unbalanced data (missing visits) naturally

# ===============================================
# 2: Fitting a basic model
# ===============================================

# The mmrm formula has two parts:
# 1. Fixed effects: FEV1 ~ RACE + SEX + ARMCD * AVISIT
# 2. Covariance specification: us(AVISIT | USUBJID) *ask more about this
#   - mean that it models the within-subject correlation using an unstructured cov matrix
#     where AVISIT defiens the repeated time points and USUBJID identifies the subjects

fit = mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)

# Check object class
class(fit)
# [1] "mmrm"     "mmrm_fit" "mmrm_tmb"
# mmrm is the top level user facing model
# mmrm_fit is a mid level class
# mmrm_tmb shows the model was fit using Template Model Builder (TMB)

# Whats inside the fitted model
names(fit)
# Important:
# 1. beta_est (coefficients) / coef(fit)
#   - 11 values
#   - effect of each predictor on FEV1
# 2. beta_vcov (variance-covariance of coefficients) / vcov(fit)
#   - 11x11 matrix (121 entries), clubSandwich?
#   - how uncertain we are about those effects
# 3. cov (the working covariance matrix) / VarCorr(fit)
#   - 4x4 matrix which shows the covariance of a subject's repeated measurements
#   - how a subject's measurements correlate across different visits

summary(fit)

# ===============================================
# 3: Extracting the pieces of the model
# ===============================================

# Fixed effects coefficients
coef(fit)

# Variance-covariance of the coefficients (not the working covariance)
vcov(fit)
# off-diagonal values = 0 means the two coefficients are independent

diag(vcov(fit))
# diag shows the variances of each coefficient - bigger diagonal value = more uncertainty for the coefficient

# Design matrix
X = model.matrix(fit)
dim(X)  # Has fewer rows than the original fev_data as it drops the rows with NA
# Each row is one observation, each column is one predictor, and the intercept column is all 1's

# Residuals
r = residuals(fit)
# observed FEV1 - predicted FEV1
length(r)  # Also fewer rows than fev_data as it drops rows with NA

# Fitted values
head(fitted(fit))
# Basically the predicted values

# Model frame: the data actually used in fitting
mf = model.frame(fit)
dim(mf)  # Again, any row with NA is dropped. 
# Also, only the variables that appear in the formula are kep 

# Formula
formula(fit)

# nobs() doesnt work yet, should return number of observations
tryCatch(nobs(fit), error = function(e) message("nobs failed: ", e$message))
# But it can be extracted simply by doing component(fit, "n_obs")

# Checking the weights
weights(fit)  # All are 1's meaning every observations contributes equally to the model

# ===============================================
# 4: Working Covariance (main part of the mmrm)
# ===============================================

# VarCorr extracts the estimated working covariance matrix
vc = VarCorr(fit)
# 4x4 matrix
dim(vc)


# To get the covariance for a subject with only 3 out of 4 visits,
# simple drop the row and columns for the missing visit

# Example: subject "PT1", which visits do they have?
subset(fev_data, USUBJID == "PT1")[, c("USUBJID", "AVISIT", "FEV1")]
# If PT1 has visits 1,2,4 but not 3, their working covariance would be:
# vc[c(1,2,4), c(1,2,4)]

fit$cov  # same as VarCorr(fit)

# ===============================================
# 5: Parts of formula: how mmrm parses the covariance spec
# ===============================================

fit$formula_parts

fit$formula_parts$formula # basically just the original formula

fit$formula_parts$model_formula # fixed effects only formula
# ARMCD:AVISIT -> ARMCD + AVISIT + ARMCD:AVISIT

fit$formula_parts$full_formula # same as $model_formula, just added USUBJID

fit$formula_parts$cov_type # the covariance structure which in this case is unstructured ("us")

fit$formula_parts$is_spatial

fit$formula_parts$visit_var  # which column identifies the time points
# Extracted from the left side of AVISIT | USUBJID

fit$formula_parts$subject_var # which column identifies the subjects
# Extracted from the right side of AVISIT | USUBJID

fit$formula_parts$model_var # predictor variables used in the fixed effects part of the model

fit$formula_parts$response_var # outcome variable


# ===============================================
# 6: component() accessor
# ===============================================

# mmrm official way to extract internals
component(fit, "n_obs") # number of total observations
component(fit, "n_subjects") # number of unique subjects
component(fit, "n_timepoints") # number of unique visits
component(fit, "subject_var") # same as fit$formula_parts$subject_va

# The "full_frame" is the data used in fitting with NAs removed
ff = component(fit, "full_frame")
dim(ff)
head(ff)
names(ff) # includes USUBJID

# The TMB level data has even more detail
names(fit$tmb_data)
fit$tmb_data$n_visits # same as component(fit, "n_timepoints")
fit$tmb_data$n_subjects # same as component(fit, "n_subjects")
head(fit$tmb_data$subject_n_visits) # visits per subject in the fitted data
head(fit$tmb_data$coordinates) # 0-based visit indices per observation, 0 = visit 1


# ===============================================
# 7: Different covariance structures
# ===============================================

fit_us = mmrm(FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID), data = fev_data) # unstructured
fit_cs = mmrm(FEV1 ~ ARMCD * AVISIT + cs(AVISIT | USUBJID), data = fev_data) # compound symmetry, one variance + one covariance
fit_ar1 = mmrm(FEV1 ~ ARMCD * AVISIT + ar1(AVISIT | USUBJID), data = fev_data) # correlation decays with distance
fit_toep = mmrm(FEV1 ~ ARMCD * AVISIT + toep(AVISIT | USUBJID), data = fev_data) # Toeplitz

VarCorr(fit_us) # all entries free, no restriction
VarCorr(fit_cs) # constant diagonal and off diagonal
VarCorr(fit_ar1) # correlation decays with lag
VarCorr(fit_toep) # each lag gets its own value

# Number of parameters, more parameters = more flexible
length(fit_us$theta_est)
length(fit_cs$theta_est)
length(fit_ar1$theta_est)
length(fit_toep$theta_est)


# ===============================================
# 8: Weighted models
# ===============================================

# mmrm supports observation level weights
fev_data$w = runif(nrow(fev_data), 0.5, 2) # just random fake weight numbers

fit_w = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  weights = fev_data$w
)

# weights are no longer all 1's
weights(fit_w)[1:10]

# The coefficients and VarCorr are different from the unweighted model
cbind(unweighted = coef(fit_us), weighted = coef(fit_w))
VarCorr(fit_w)
VarCorr(fit_us)


# ===============================================
# 9: Degrees of freedom methods
# ===============================================

fit_satt = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  method = "Satterthwaite"
)

fit_kr = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  method = "Kenward-Roger"
)

# Compare the coefficient tables
summary(fit_satt)$coefficients
summary(fit_kr)$coefficients

# will need to look into more of both methods and other methods mentioned in the website


