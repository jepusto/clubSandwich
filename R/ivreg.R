#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an ivreg object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an ivreg object fitted
#' from the \CRANpkg{AER} package or the \CRANpkg{ivreg} package.
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
#' @details For any "ivreg" objects fitted via the \code{\link[ivreg]{ivreg}} 
#'   function from the \CRANpkg{ivreg} package, only traditional 2SLS 
#'   regression method (method = "OLS") is supported.
#'   clubSandwich currently cannot support robust-regression methods such as
#'   M-estimation (method = "M") or MM-estimation (method = "MM").
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists
#'   of a matrix of the estimated variance of and covariances between the
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'
#' @examples 
#' 
#' if (requireNamespace("AER", quietly = TRUE)) withAutoprint({
#' 
#'   library(AER)
#'   data("CigarettesSW")
#'   Cigs <- within(CigarettesSW, {
#'     rprice <- price/cpi
#'     rincome <- income/population/cpi
#'     tdiff <- (taxs - tax)/cpi
#'   })
#' 
#'   iv_fit_AER <- AER::ivreg(log(packs) ~ log(rprice) + log(rincome) | 
#'                   log(rincome) + tdiff + I(tax/cpi), data = Cigs)
#'   vcovCR(iv_fit_AER, cluster = Cigs$state, type = "CR2")
#'   coef_test(iv_fit_AER, vcov = "CR2", cluster = Cigs$state)
#'
#' })
#' 
#' pkgs_available <- 
#'   requireNamespace("AER", quietly = TRUE) & 
#'   requireNamespace("ivreg", quietly = TRUE)
#'
#' if (pkgs_available) withAutoprint ({
#' 
#' data("CigarettesSW")
#'   Cigs <- within(CigarettesSW, {
#'     rprice <- price/cpi
#'     rincome <- income/population/cpi
#'     tdiff <- (taxs - tax)/cpi
#'   })
#' iv_fit_ivreg <- ivreg::ivreg(log(packs) ~ log(rprice) + log(rincome) | 
#'                   log(rincome) + tdiff + I(tax/cpi), data = Cigs)
#'   vcovCR(iv_fit_ivreg, cluster = Cigs$state, type = "CR2")
#'   coef_test(iv_fit_ivreg, vcov = "CR2", cluster = Cigs$state)
#' })
#' 
#' @export

vcovCR.ivreg <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (inverse_var != FALSE) stop("Unfortunately, the inverse_var option is not available for ivreg models.")
  if (!is.null(obj$method) && obj$method %in% c("M", "MM")) stop("clubSandwich does not currently support ivreg models estimated using method = 'M' or method = 'MM'.")
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

#' @export

model_matrix.ivreg <- function(obj) {
  model_matrix <- model.matrix(obj, component = "projected")
  
  w <- obj$weights
  if (is.null(w) || all(pos_wts <- w > 0)) {
    return(model_matrix)
  } else {
    return(model_matrix[pos_wts > 0,,drop=FALSE])
  }
  
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.ivreg() is in AER package
# use default v_scale()