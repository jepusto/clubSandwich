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
