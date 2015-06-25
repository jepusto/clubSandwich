# residuals_CR()
# coef()
# nobs()

#-------------------------------------
# model_matrix()
#-------------------------------------

model_matrix.gls <- function(obj) {
  model.matrix(formula(obj), data = nlme::getData(obj))
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.gls <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  V <- matrix(0, N, N)
  for (i in levels(groups)) {
    V[groups == i, groups == i] <- nlme::getVarCov(obj, individual = i)
  }
  V
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.gls <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  W <- matrix(0, N, N)
  for (i in levels(groups)) {
    W[groups == i, groups == i] <- chol2inv(chol(nlme::getVarCov(obj, individual = i)))
  }
  W
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.gls <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- nlme::getGroups(obj)
  if (missing(target)) target <- NULL
  if (missing(inverse_var) ) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}