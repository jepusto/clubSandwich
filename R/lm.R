# model_matrix()
# residuals_CR()
# coef()

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.lm <- function(obj) {
  rep(1, nobs(obj))
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.lm <- function(obj) {
  weights <- weights(obj)
  if (is.null(weights)) weights <- 1
  rep(weights, length.out = nobs(obj))
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.lm <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(weights) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}
