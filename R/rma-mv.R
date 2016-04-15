
#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a robu object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a 
#' \code{\link[metafor]{rma.mv}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to the factor in
#'   the random-effects structure with the fewest distinct levels. Caveat
#'   emptor: the function does not check that the random effects are nested.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the 
#'   estimated variance-covariance structure of the \code{rma.mv} object.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export
#' 
#' @examples
#' library(metafor)
#' data(hierdat, package = "robumeta")
#' 
#' mfor_fit <- rma.mv(effectsize ~ binge + followup + sreport + age, 
#'                  V = var, random = list(~ 1 | esid, ~ 1 | studyid),
#'                  data = hierdat)
#' mfor_fit
#' 
#' mfor_CR2 <- vcovCR(mfor_fit, type = "CR2")
#' mfor_CR2
#' coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
#' 
#' Wald_test(mfor_fit, constraints = c(2,4), vcov = mfor_CR2)
#' Wald_test(mfor_fit, constraints = 2:5, vcov = mfor_CR2)

vcovCR.rma.mv <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- findCluster.rma.mv(obj)
  if (missing(target)) {
    target <- NULL
    inverse_var <- is.null(obj$W)
  } else {
    if (missing(inverse_var)) inverse_var <- FALSE
  }
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

# coef()
# residuals_CR()
# vcov()
# model_matrix

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.rma.mv <- function(obj, cluster) {
  matrix_list(obj$M, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.rma.mv <- function(obj, cluster) {
  if (is.null(obj$W)) {
    V_list <- targetVariance(obj, cluster)
    lapply(V_list, function(v) chol2inv(chol(v)))
  } else{
    matrix_list(obj$W, cluster, "both")
  }
}

#-----------------------------------------------
# Get outer-most clustering variable
#-----------------------------------------------

findCluster.rma.mv <- function(obj) {
  if (obj$withS) {
    r <- which.min(obj$s.nlevels)
    cluster <- obj$mf.r[[r]][[obj$s.names[r]]]
  } else if (obj$withG) {
    cluster <- obj$mf.r[[1]][[obj$g.names[2]]]
  } else {
    stop("No clustering variable specified.")
  }
  droplevels(as.factor(cluster))
}
