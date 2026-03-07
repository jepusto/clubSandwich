context("mmrm objects")
set.seed(20190513)

skip_if_not_installed("mmrm")
skip_if_not_installed("nlme")

library(mmrm, quietly = TRUE, warn.conflicts = FALSE)

# Dataset 1: fev_data (from mmrm package)

data(fev_data, package = "mmrm")

fev_data$wt <- 1L + rpois(nrow(fev_data), lambda = 3)

# Unstructured covariance
obj_us <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
  data = fev_data,
  weights = fev_data$wt
)

# Toeplitz covariance
obj_toep <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + toep(AVISIT | USUBJID),
  data = fev_data
)

# Grouped covariance (separate covariance per treatment arm)
obj_grouped <- mmrm(
  FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | ARMCD / USUBJID),
  data = fev_data
)

# Extract cluster variables for fev_data models
ff_us <- component(obj_us, "full_frame")
cluster_us <- droplevels(as.factor(ff_us[[obj_us$formula_parts$subject_var]]))

ff_toep <- component(obj_toep, "full_frame")
cluster_toep <- droplevels(as.factor(ff_toep[[obj_toep$formula_parts$subject_var]]))

ff_grouped <- component(obj_grouped, "full_frame")
cluster_grouped <- droplevels(as.factor(ff_grouped[[obj_grouped$formula_parts$subject_var]]))


# Dataset 2: BodyWeight (from nlme) — rat growth data
# 16 rats, 3 diet groups, 11 time points

data(BodyWeight, package = "nlme")
bw_data <- as.data.frame(BodyWeight)
bw_data$Time_f <- factor(bw_data$Time)

# Compound symmetry covariance
obj_bw <- mmrm(
  weight ~ Diet + cs(Time_f | Rat),
  data = bw_data
)

ff_bw <- component(obj_bw, "full_frame")
cluster_bw <- droplevels(as.factor(ff_bw[[obj_bw$formula_parts$subject_var]]))


# Dataset 3: Orthodont (from nlme) — orthodontic growth data
# 27 subjects (16 Male, 11 Female), ages 8/10/12/14

data(Orthodont, package = "nlme")
orth_data <- as.data.frame(Orthodont)
orth_data$age_f <- factor(orth_data$age)
orth_data$Subject <- factor(as.character(orth_data$Subject))

# AR(1) covariance (interaction breaks perfect orthogonality in design)
obj_orth <- mmrm(
  distance ~ Sex * age_f + ar1(age_f | Subject),
  data = orth_data
)

ff_orth <- component(obj_orth, "full_frame")
cluster_orth <- droplevels(as.factor(ff_orth[[obj_orth$formula_parts$subject_var]]))
