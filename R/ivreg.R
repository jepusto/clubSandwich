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
#' @param inverse_var Not used for \code{ivreg} objects.
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
#' data("CigarettesSW", package = "AER")
#' Cigs <- within(CigarettesSW, {
#'   rprice <- price/cpi
#'   rincome <- income/population/cpi
#'   tdiff <- (taxs - tax)/cpi
#' })
#' 
#' iv_fit <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi), data = Cigs)
#' vcovCR(iv_fit, cluster = Cigs$state, type = "CR2")
#' coef_test(iv_fit, vcov = "CR2", cluster = Cigs$state)
#'       
#' @export

vcovCR.ivreg <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (inverse_var != FALSE) stop("Unfortunately, the inverse_var option is not available for ivreg models.")
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# residuals_CS()
# coef()
# targetVariance()
# weightMatrix()
# v_scale()

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

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.ivreg() is in AER package
# use default v_scale()