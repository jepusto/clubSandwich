#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an lm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{\link{lm}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{lm} objects.
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

vcovCR.lm <- function(obj, cluster, type, target = NULL, inverse_var = NULL) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

# model_matrix()
# residuals_CS()
# coef()
# nobs()
# targetVariance()
# weightMatrix()


#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.lm() is in sandwich package

v_scale.lm <- function(obj) {
  as.vector(sum(summary(obj)$df[1:2]))
}