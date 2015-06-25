# nobs()

#-------------------------------------
# residuals_CR()
#-------------------------------------

residuals_CR.lme <- function(obj) 
  residuals(obj, level = 0)

#-------------------------------------
# coef_CR()
#-------------------------------------

coef_CR.lme <- function(obj)
  nlme::fixef(obj)

#-------------------------------------
# model_matrix()
#-------------------------------------

model_matrix.lme <- function(obj) {
  model.matrix(formula(obj), data = nlme::getData(obj))
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.lme <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  V <- matrix(0, N, N)
  for (i in levels(groups)) {
    V[groups == i, groups == i] <- nlme::getVarCov(obj, individual = i, type = "marginal")[[1]]
  }
  V
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.lme <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  W <- matrix(0, N, N)
  for (i in levels(groups)) {
    W[groups == i, groups == i] <- chol2inv(chol(nlme::getVarCov(obj, individual = i, type = "marginal")[[1]]))
  }
  W
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.lme <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- nlme::getGroups(obj)
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}