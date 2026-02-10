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
# 2. Covariance spec: us(AVISIT | USUBJID)

fit = mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
  data = fev_data
)

# Check object class
class(fit)
# [1] "mmrm"     "mmrm_fit" "mmrm_tmb"

# Whats inside the fitted model
names(fit)

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
