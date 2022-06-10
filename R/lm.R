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
#' @examples 
#' 
#' data("ChickWeight", package = "datasets")
#' lm_fit <- lm(weight ~ Time + Diet:Time, data = ChickWeight)
#' vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#' 
#' if (requireNamespace("plm", quietly = TRUE)) {
#' 
#'   data("Produc", package = "plm")
#'   lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#'   individual_index <- !grepl("state", names(coef(lm_individual)))
#'   vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index]
#' 
#'   # compare to plm()
#'   plm_FE <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#'                      data = Produc, index = c("state","year"), 
#'                      effect = "individual", model = "within")
#'   vcovCR(plm_FE, type="CR2")
#'   
#' }
#' 
#' @export

vcovCR.lm <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
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

#' @export

v_scale.lm <- function(obj) {
  as.vector(sum(summary(obj)$df[1:2])) 
}
