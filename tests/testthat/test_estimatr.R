context("estimatr objects")
set.seed(20190513)

skip_if_not_installed("estimatr")

library(estimatr)

data(mtcars)

# =============== vcovCR ===============

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)

vcov_lm <- vcovCR(lm_fit, ChickWeight$Chick, "CR2") # works
vcov_lmr <- vcovCR(lm_robust, ChickWeight$Chick, "CR2") # doesn't work

# expect_equal(vcov_lm, vcov_lmr)

# =============== model_matrix() ===============

mm_fit <- model_matrix(lm_fit) # works
mm_rob <- model_matrix(lm_rob) # doesn't work

# expect_equal(mm_fit, mm_rob)


# =============== residuals_CS() ===============

rcs_fit <- residuals_CS(lm_fit) # works
rcs_rob <- residuals_CS(lm_rob) # gives null

# expect_equal(rcs_fit, rcs_rob)

# =============== coef() ===============

coef_fit <- coef(lm_fit) # works
coef_rob <- coef(lm_rob) # works?

expect_equal(coef_fit, coef_rob) # works!!

# =============== nobs() ===============

nobs_fit <- nobs(lm_fit) # works
nobs_rob <- nobs(lm_rob) # works

expect_equal(nobs_fit, nobs_rob) # works!!!

# =============== targetVariance() ===============

tV_fit <- targetVariance(lm_fit, ChickWeight$Chick) # works
tV_rob <- targetVariance(lm_rob, ChickWeight$Chick) # also works

expect_equal(tV_fit, tV_rob) # works!

# =============== weightMatrix() ===============

wM_fit <- weightMatrix(lm_fit, ChickWeight$Chick) # works!!
wM_rob <- weightMatrix(lm_rob, ChickWeight$Chick) # works!

expect_equal(wM_fit, wM_rob) # works!!!!

# =============== v_scale() ===============

vs_fit <- v_scale(lm_fit) # works
vs_rob <- v_scale(lm_rob) # works

expect_equal(vs_fit, vs_rob) # works!!


# lm1 <- lm(wt ~ mpg, data = mtcars)
# lmr1 <- lm_robust(wt ~ mpg, data = mtcars)



# data(Duncan, package = "carData")
# Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#   lm2 <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
# lmr2 <- 
#   lm_robust(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
# 
# Wald1 <- Wald_test(lm2, 
#           constraints = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE), 
#           vcov = "CR2", cluster = Duncan$cluster)
# Wald2 <- Wald_test(lmr2, 
#           constraints = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE), 
#           vcov = "CR2", cluster = Duncan$cluster)
# # causes error ^
# # expect_equal(Wald1, Wald2)
# 
# 
# data("ChickWeight", package = "datasets")
# lm3 <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
# lmr3 <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
# 
# 
# 
# model <- stats::lm(mpg ~ wt + hp, data = mtcars)
# summary(model)
# 
# robust_model <- lm_robust(mpg ~ wt + hp, data = mtcars, se_type = "HC2")
# summary(robust_model)
# 
# # aer_ivereg
# # none?
# 
# # coef
# balanced_dat <- function(m, n) {
#   cluster <- factor(rep(LETTERS[1:m], each = n))
#   N <- length(cluster)
#   m1 <- sample(3:(m-7), size = 1)
#   m2 <- sample((m1 + 3):(m-3), size = 1) - m1
#   m3 <- m - m1 - m2
#   c(m1, m2, m3)
#   X_btw <- rep(rep(LETTERS[1:3], c(m1, m2, m3)), each = n)
#   X_wth <- rep(rep(c(0,1), each = n / 2), m)
#   nu <- rnorm(m)[cluster]
#   e <- rnorm(n * m)
#   y <- nu + e
#   data.frame(y, X_btw, X_wth, cluster, row = 1:N)
# }
# 
# dat <- balanced_dat(m = 15, n = 8)
# lm_fit <- lm(y ~ X_btw + X_wth, data = dat)
# which_grid <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(lm_fit))))
# tests_all <- coef_test(lm_fit, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
# 
# tests_A <- apply(which_grid[-1,], 1, function(x) tests_all[x,])
# 
# lm_rob <- lm_robust(y ~ X_btw + X_wth, data = dat)
# which_gridB <- expand.grid(rep(list(c(FALSE,TRUE)), length(coef(lm_rob))))
# tests_allB <- coef_test(lm_rob, vcov = "CR0", cluster = dat$cluster, test = "All", coefs = "All", p_values = FALSE)
# 
# tests_B <- apply(which_gridB[-1,], 1, function(x) tests_allB[x,])
# 
# # conf_int
# data("ChickWeight", package = "datasets")
# lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
# lm_rob <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
# lc1 <- linear_contrast(lm_fit,
#                        vcov = "CR2",
#                        cluster = ChickWeight$Chick, 
#                        contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
#                        p_values = TRUE)
# lc2 <- linear_contrast(lm_rob,
#                        vcov = "CR2",
#                        cluster = ChickWeight$Chick, 
#                        contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE),
#                        p_values = TRUE)
# lm_rob2 <- lm_robust(weight ~ 0 + Diet + Time:Diet, data = ChickWeight,
#                      clusters = Chick)

# estfun


# geeglm


# glm_logit