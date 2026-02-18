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

# Variance-covariance of the coefficients 
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
# 4: Working covariance (main part of the mmrm)
# ===============================================

# VarCorr extracts the estimated working covariance matrix
vc = VarCorr(fit)
# 4x4 matrix
dim(vc)


# To get the covariance for a subject with only 3 out of 4 visits,
# simply drop the row and columns for the missing visit

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


# ===============================================
# mmrm 'sandwich' estimators
# ===============================================

# mmrm has built in empirical (sandwich) variance estimators

fit_cr0 = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  control = mmrm_control(vcov = "Empirical")
)
# meat = Σ_j X_j 'e_j e_j' X_j
# problem: biased estimate of true residual covariance as it uses B_hat and not true B

fit_cr2 = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  control = mmrm_control(vcov = "Empirical-Bias-Reduced")
)
# meat = Σ_j X_j 'A_j e_j e_j' A_j' X_j (adds the A_j correction matrix)
# approximately unbiased

fit_cr3 = mmrm(
  FEV1 ~ ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data,
  control = mmrm_control(vcov = "Empirical-Jackknife")
)
# meat = Σ_j X_j' (I - H_jj)^(-1) e_j e_j' (I - H_jj)^(-1) X_j (adds H_jj, hat matrix)



# Compare coefficient standard errors across methods
se_asymptotic = summary(fit_us)$coefficients[, "Std. Error"]
se_cr0 = summary(fit_cr0)$coefficients[, "Std. Error"]
se_cr2 = summary(fit_cr2)$coefficients[, "Std. Error"]
se_cr3 = summary(fit_cr3)$coefficients[, "Std. Error"]

data.frame(
  asymptotic = se_asymptotic,
  CR0 = se_cr0,
  CR2 = se_cr2,
  CR3 = se_cr3
)
# CR0 < CR2 < CR3


# ===============================================
# Potential clubSandwich integration
# ===============================================

# clubSandwich needs to construct these for each subject (cluster):
#   1. X_j: design matrix rows for subject j
#   2. W_j: weight matrix for subject j (inverse of working covariance)
#   3. Theta_j: target variance (working covariance block for subject j)
#   4. e_j: residuals for subject j


# Manually building them 
ff = component(fit_us, "full_frame")
subj_var = fit_us$formula_parts$subject_var
visit_var = fit_us$formula_parts$visit_var
visit_levels = levels(ff[[visit_var]])

subjects = ff[[subj_var]]
visits = ff[[visit_var]]
vc_full = VarCorr(fit_us)

# For one subject, build the working covariance sub-matrix:
subj1 = levels(subjects)[1]
subj1_visits = visits[subjects == subj1]
subj1_idx = match(as.character(subj1_visits), visit_levels)

cat("Subject:", subj1, "\n")
cat("Visits:", as.character(subj1_visits), "\n")
cat("Visit indices:", subj1_idx, "\n")
cat("Working covariance sub-matrix:\n")
print(vc_full[subj1_idx, subj1_idx])

# To get Theta_j: target variance (working covariance block for subject j)
# Split observation indices by subject
subjects = droplevels(ff[[subj_var]])
obs_by_subject = split(seq_along(subjects), subjects)

# For each subject, get their visit indices and extract the sub-matrix
target_list = lapply(obs_by_subject, function(rows) {
  subj_visits = visits[rows]
  idx = match(as.character(subj_visits), visit_levels)
  vc_full[idx, idx]
})
# loop over each subject's rows, looks up their visits, converts to indices, and extracts the sub matrix

target_list[11]


# To get X_j — design matrix rows per subject:
X = model.matrix(fit_us)

# For PT1 
X[obs_by_subject[["PT1"]], ] # model_matrix.mmrm()
# A 2×11 matrix: PT1's 2 observations, 8 coefficients


# To get W_j - weight matrix per subject (inverse of working covariance):
weight_list = lapply(target_list, solve) # solve() computes the matrix inverse

# For PT1:
weight_list[[1]] # weightMatrix.mmrm() 
# A 3×3 matrix — the inverse of PT1's covariance block


# To get e_j — residuals per subject:
e = residuals(fit_us)

# For PT1:
e[obs_by_subject[["PT1"]]] # residuals_CS.mmrm()


# function that takes a fitted mmrm object and returns a list of per-subject working covariance matrices.

my_targetVariance = function(obj) {
  ff = component(obj, "full_frame")
  subj_var = obj$formula_parts$subject_var
  visit_var = obj$formula_parts$visit_var
  visit_levels = levels(ff[[visit_var]])
  
  subjects = ff[[subj_var]]
  visits = ff[[visit_var]]
  vc_full = VarCorr(obj)
  
  obs_by_subject = split(seq_along(subjects), subjects)
  
  lapply(obs_by_subject, function(rows) {
    idx = match(as.character(visits[rows]), visit_levels)
    vc_full[idx, idx]
  })
}


# ===============================================
# Group specific covariance
# ===============================================

# mmrm can fit separate covariance matrices per group
# Syntax: us(AVISIT | ARMCD / USUBJID)

fit_grouped = mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | ARMCD / USUBJID),
  data = fev_data
)

VarCorr(fit_grouped)

class(VarCorr(fit_grouped)) # list of matrices

fit_grouped$formula_parts$group_var # "ARMCD"
component(fit_grouped, "n_groups") # 2

# For non group case, building a subjects covariance block was:
# 
# 1. Get the full Σ (one matrix)
# 2. Find the subject's visit indices
# 3. Subset: Σ[idx, idx]
# 
# For grouped case:
# 
# 1. Get the list of Σs (one per group)
# 2. Figure out which group this subject belongs to
# 3. Pick the correct Σ from the list
# 4. Find the subject's visit indices
# 5. Subset: Σ_group[idx, idx]


# Figuring out each subject's group:
ff = component(fit_grouped, "full_frame")
group_var = fit_grouped$formula_parts$group_var # "ARMCD"
groups = ff[[group_var]] # factor: (PBO, PBO, PBO, TRT, TRT, etc)
vc_list = VarCorr(fit_grouped)  # list: $PBO = 4×4, $TRT = 4×4

# obs_by_subject = split(seq_along(subjects), droplevels(subjects))
#
# lapply(obs_by_subject, function(rows) {
#   subj_group  as.character(groups[rows[1]]) # "PBO" or "TRT"
#   vc_this  vc_list[[subj_group]] # pick the right Σ
#   idx  match(as.character(visits[rows]), visit_levels)
#   vc_this[idx, idx] # subset by visits
# })

# Updated function would be:
targetVariance.mmrm = function(obj, cluster) {
  ff = component(obj, "full_frame")
  visit_var = obj$formula_parts$visit_var
  group_var = obj$formula_parts$group_var
  visit_levels = levels(ff[[visit_var]])
  visits = ff[[visit_var]]
  vc = VarCorr(obj)
  
  obs_by_subject = split(seq_along(cluster), cluster)
  
  if (is.null(group_var)) {
    # Non grouped: one Σ for everyone
    lapply(obs_by_subject, function(rows) {
      idx = match(as.character(visits[rows]), visit_levels)
      vc[idx, idx]
    })
  } else {
    # Grouped: pick the right Σ per subject
    groups = ff[[group_var]]
    lapply(obs_by_subject, function(rows) {
      subj_group = as.character(groups[rows[1]])
      idx = match(as.character(visits[rows]), visit_levels)
      vc[[subj_group]][idx, idx]
    })
  }
}
