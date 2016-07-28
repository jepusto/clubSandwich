#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an ivreg object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{\link[AER]{ivreg}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{ivreg} objects.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be an identity matrix.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists
#'   of a matrix of the estimated variance of and covariances between the
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

vcovCR.ivreg <- function(obj, cluster, type, target = NULL, inverse_var = NULL) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

# residuals_CS()
# coef()

#----------------------------------------------
# get X matrix
#----------------------------------------------

model_matrix.ivreg <- function(obj) {
  model.matrix(obj, component = "regressors")
}

#----------------------------------------------
# get projection matrix
#----------------------------------------------

projection_matrix.ivreg <- function(obj) {
  model.matrix(obj, component = "projected")
}

#--------------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------------

targetVariance.ivreg <- function(obj) {
  rep(1, nobs(obj))
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.ivreg <- function(obj) {
  weights <- weights(obj)
  if (is.null(weights)) weights <- 1
  rep(weights, length.out = nobs(obj))
}
