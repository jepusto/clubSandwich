#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an lme object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a \code{\link[nlme]{lme}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to 
#'   \code{getGroups(obj)}.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the
#'   estimated variance-covariance structure of the \code{lme} object.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

vcovCR.lme <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- nlme::getGroups(obj)
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

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
