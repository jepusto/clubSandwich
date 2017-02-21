#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a robu object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a
#' \code{\link[robumeta]{robu}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to the
#'   \code{studynum} used in fitting the \code{\link[robumeta]{robu}} object.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the 
#'   inverse of the estimated weights used in fitting the
#'   \code{\link[robumeta]{robu}} object.
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
#' library(robumeta)
#' data(hierdat)
#' 
#' robu_fit <- robu(effectsize ~ binge + followup + sreport + age, 
#'                  data = hierdat, studynum = studyid, 
#'                  var.eff.size = var, modelweights = "HIER")
#' robu_fit
#' 
#' robu_CR2 <- vcovCR(robu_fit, type = "CR2")
#' robu_CR2
#' coef_test(robu_fit, vcov = robu_CR2, test = c("Satterthwaite", "saddlepoint"))
#' 
#' Wald_test(robu_fit, constraints = c(2,4), vcov = robu_CR2)
#' Wald_test(robu_fit, constraints = 2:5, vcov = robu_CR2)
#' 



vcovCR.robu <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) cluster <- obj$study_orig_id
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target) & (!obj$user_weighting)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

#-----------------------------------------------
# coefficients
#-----------------------------------------------

coef_CS.robu <- function(obj) {
  beta <- as.vector(obj$b.r)
  labs <- obj$reg_table$labels
  names(beta) <- levels(labs)[labs]
  beta
}

#-----------------------------------------------
# residuals
#-----------------------------------------------

residuals_CS.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  obj$data.full$e.r[ord]
}


#-----------------------------------------------
# Model matrix
#-----------------------------------------------

model_matrix.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  obj$Xreg[ord,]
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.robu <- function(obj, cluster) {
  ord <- order(order(obj$study_orig_id))
  if (obj$user_weighting) {
    V <- obj$data.full$avg.var.eff.size[ord]
  } else {
    V <- 1 / obj$data.full$r.weights[ord]
  }
  matrix_list(V, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.robu <- function(obj, cluster) {
  ord <- order(order(obj$study_orig_id))
  if (obj$user_weighting) { 
    W <- obj$data.full$userweights[ord]
  } else{
    W <- obj$data.full$r.weights[ord]
  }
  W <- W # / mean(wi)
  matrix_list(W, cluster, "both")
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

bread.robu <- function(x, ...) {
  if (x$user_weighting) { 
    W <- x$data.full$userweights
  } else{
    W <- x$data.full$r.weights
  }
  x$N * chol2inv(chol(crossprod(x$Xreg, W * x$Xreg)))
}

v_scale.robu <- function(obj) {
  obj$N
}
