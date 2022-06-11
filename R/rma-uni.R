#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a rma.uni object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a 
#' \code{\link[metafor]{rma.uni}} object.
#' 
#' @param cluster Expression or vector indicating which observations 
#'   belong to the same cluster. Required for \code{rma.uni} objects.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be diagonal
#'   with entries equal to the estimated marginal variance of the effect sizes. 
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
#' if (requireNamespace("metafor", quietly = TRUE) & requireNamespace("metadat", quietly = TRUE)) {
#' 
#' library(metafor)
#' data(dat.assink2016, package = "metadat")
#' 
#' mfor_fit <- rma.uni(yi ~ year + deltype, vi = vi, 
#'                  data = dat.assink2016)
#' mfor_fit
#' 
#' mfor_CR2 <- vcovCR(mfor_fit, type = "CR2", cluster = dat.assink2016$study)
#' mfor_CR2
#' coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
#' Wald_test(mfor_fit, constraints = constrain_zero(2:4), vcov = mfor_CR2)
#' 
#' }
#' 


vcovCR.rma.uni <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {

  if (missing(cluster)) stop("You must specify a clustering variable.")
  
  # if (length(cluster) != nrow(model_matrix(obj))) cluster <- droplevels(as.factor(cluster[obj$not.na]))
  # if (length(cluster) != nrow(model_matrix(obj))) stop("Clustering variable must have length equal to nrow(model_matrix(obj)).")

  if (missing(target)) {
    target <- NULL
    if (missing(inverse_var)) inverse_var <- is.null(obj$weights) & obj$weighted
  } else {
    if (missing(inverse_var)) inverse_var <- FALSE
  }
  
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# residuals_CS()

#' @export

residuals_CS.rma <- function(obj) {
  res <- residuals(obj)
  not_na <- obj$not.na
  
  if (length(res) == length(not_na)) res <- res[not_na]
  if (!is.null(wts <- weights(obj)) && !all(pos_wts <- wts > 0)) {
    res <- res[pos_wts]
  }
  
  return(res)
}

# coef()
# vcov()
# model_matrix()

# na.action

#' @export

na.action.rma <- function(object, ...) {
  
  if (all(object$not.na)) return(NULL)
  
  res <- which(!object$not.na)
  class(res) <- "omit"
  res
  
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

#' @export

targetVariance.rma.uni <- function(obj, cluster) {
  
  vi <- obj$vi
  
  if (obj$weighted && !is.null(wts <- obj$weights) && !all(pos_wts <- wts > 0)) {
    vi <- vi[pos_wts]
  }
  
  matrix_list(vi + obj$tau2, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weightMatrix.rma.uni <- function(obj, cluster) {
  if (obj$weighted) {
    if (is.null(obj$weights)) {
      wi <- 1 / (obj$vi + obj$tau2)  
    } else {
      wi <- obj$weights
      wi <- wi[wi > 0]
    }
  } else {
    wi <- rep(1, obj$k)
  }
  w_scale <- mean(wi)
  wi <- wi / w_scale
  
  W_list <- matrix_list(wi, cluster, "both")
  attr(W_list, "w_scale") <- w_scale
  W_list
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.rma.uni <- function(x, ...) {
  X_mat <- model_matrix(x)
  if (x$weighted) {
    if (is.null(x$weights)) {
      wi <- 1 / (x$vi + x$tau2)  
    } else {
      wi <- x$weights
      wi <- wi[wi > 0]
    }
    XWX <- crossprod(X_mat, wi * X_mat)
  } else {
    XWX <- crossprod(X_mat)
  }
  B <- chol2inv(chol(XWX)) * nobs(x)
  rownames(B) <- colnames(B) <- colnames(X_mat)
  B
}

#' @export

v_scale.robu <- function(obj) {
  nobs(obj)
}
