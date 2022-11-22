#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a fixest object.
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates from an \code{\link{fixest}}
#' object estimated using \code{method = "feols"}.
#'
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. If omitted, it will be inferred from the structure of \code{obj}. 
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
#' @examples
#'
#'
#' if (requireNamespace("fixest", quietly = TRUE)) withAutoprint({
#'
#' data("base_did")
#' est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
#' vcovCR(est_did, type = "CR2")
#'
#'   # compare to feols() clustered SEs
#'
#' })
#'
#' @export

vcovCR.fixest <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}


#-------------------------------------
# model_matrix()
#-------------------------------------

#' @export

model_matrix.fixest <- function(obj) {
 model.matrix(obj, type = "rhs") 
}

#-------------------------------------
# coef_CS()
#-------------------------------------

# Use fixest::coef.fixest() via default method

#-------------------------------------
# residuals_CS()
#-------------------------------------

# Use fixest::residuals.fixest() via default method

#-------------------------------------
# nobs()
#-------------------------------------

#' Use fixest::nobs.fixest()

#-------------------------------------
# targetVariance()
#-------------------------------------


#-------------------------------------
# weightMatrix()
#-------------------------------------


#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.fixest() is in the fixest package

#' @export

v_scale.fixest <- function(obj) {
  as.vector(sum(summary(obj)$df[1:2])) 
}
