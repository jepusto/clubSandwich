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
#' 
#' if (requireNamespace("robumeta", quietly = TRUE)) {
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
#' Wald_test(robu_fit, constraints = constrain_zero(c(2,4)), vcov = robu_CR2)
#' Wald_test(robu_fit, constraints = constrain_zero(2:5), vcov = robu_CR2)
#' }
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

#' @export

coef_CS.robu <- function(obj) {
  beta <- as.vector(obj$b.r)
  labs <- obj$reg_table$labels
  if (is.factor(labs)) labs <- levels(labs)[labs]
  names(beta) <- labs
  beta
}

#-----------------------------------------------
# residuals
#-----------------------------------------------

#' @export

residuals_CS.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  resid <- obj$data.full$e.r[ord]
    
  if (obj$user_weighting) {
    pos_wts <- obj$data.full$userweights[ord] > 0
    if (!all(pos_wts)) resid <- resid[pos_wts]
  } 
  
  return(resid)
}


#-----------------------------------------------
# Model matrix
#-----------------------------------------------

#' @export

model_matrix.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  model_matrix <- obj$Xreg[ord,,drop=FALSE]
  
  if (obj$user_weighting) {
    pos_wts <- obj$data.full$userweights[ord] > 0
    if (!all(pos_wts)) model_matrix <- model_matrix[pos_wts,,drop=FALSE]
  }
  
  return(model_matrix)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

#' @export

targetVariance.robu <- function(obj, cluster) {
  ord <- order(order(obj$study_orig_id))
  if (obj$user_weighting) {
    pos_wts <- obj$data.full$userweights[ord] > 0
    V <- obj$data.full$avg.var.eff.size[ord][pos_wts]
  } else {
    V <- mean(obj$data.full$r.weights) / obj$data.full$r.weights[ord]
  }
  matrix_list(V, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weights.robu <- function(object, ...) {
  ord <- order(order(object$study_orig_id))
  if (object$user_weighting) { 
    object$data.full$userweights[ord]
  } else{
    NULL
  }
}

#' @export

weightMatrix.robu <- function(obj, cluster) {
  ord <- order(order(obj$study_orig_id))
  if (obj$user_weighting) { 
    W <- obj$data.full$userweights[ord]
    W <- W[W > 0]
  } else{
    W <- obj$data.full$r.weights[ord]
  }
  w_scale <- mean(W)
  W <- W / w_scale
  W_list <- matrix_list(W, cluster, "both")
  attr(W_list, "w_scale") <- w_scale
  W_list
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.robu <- function(x, ...) {
  if (x$user_weighting) { 
    W <- x$data.full$userweights
  } else{
    W <- x$data.full$r.weights
  }
  x$N * chol2inv(chol(crossprod(x$Xreg, W * x$Xreg)))
}

#' @export

v_scale.robu <- function(obj) {
  obj$N
}
