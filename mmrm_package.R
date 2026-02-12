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
# 1. beta_est (coefficients)
# 2. beta_vcov (variance-covariance of coefficients)
# 3. cov (the working covariance matrix)

summary(fit)

# ===============================================
# 3: Extracting the pieces of the model
# ===============================================

# Fixed effects coefficients
coef(fit)

# Variance-covariance of the coefficients (NOT the working covariance)
vcov(fit)

# Design matrix
X = model.matrix(fit)
dim(X)  # Has fewer rows than the original fev_data as it drops the rows with NA

# Residuals
r = residuals(fit)
length(r)  # Also fewer rows than fev_data as it drops rows with NA

# Fitted values
head(fitted(fit))

# Model frame: the data actually used in fitting
mf = model.frame(fit)
dim(mf)  # Why is the dimension different compared to fev_data?

# Formula
formula(fit)

# nobs() doesnt work yet
tryCatch(nobs(fit), error = function(e) message("nobs failed: ", e$message))

# Checking the weights
weights(fit)  # Why are they all 1's?

# ===============================================
# 4: Working Covariance (main part of the mmrm)
# ===============================================

# VarCorr extracts the estimated working covariance matrix
vc = VarCorr(fit)
vc

dim(vc)

# Is this per-subject or for the whole dataset?
# How would you get the covariance for a subject with only 3 out of 4 visits?

# Try it: subject "PT1" — which visits do they have?
subset(fev_data, USUBJID == "PT1")[, c("USUBJID", "AVISIT", "FEV1")]

# If PT1 has visits 1,2,4 but not 3, their working covariance would be:
# vc[c(1,2,4), c(1,2,4)]

fit$cov  # same as VarCorr(fit)

# ===============================================
# 5: Parts of formula: how mmrm parses the covariance spec
# ===============================================

fit$formula_parts

# What is the subject variable? The visit variable? The covariance type?
fit$formula_parts$subject_var
fit$formula_parts$visit_var
fit$formula_parts$cov_type

# Is there a group variable? (for group-specific covariance)
fit$formula_parts$group_var


