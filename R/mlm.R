#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an mlm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{mlm} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, each row of the data will be
#'   treated as a separate cluster.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be an identity matrix.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @examples
#' iris_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + Petal.Length + Petal.Width, data = iris)
#' Vcluster <- vcovCR(iris_fit2)
#'     
#' @export

vcovCR.mlm <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  d <- ncol(residuals(obj))
  if (missing(cluster)) cluster <- 1:nobs(obj)
  if (length(cluster) == nobs(obj)) cluster <- rep(cluster, each = d)
  if (length(cluster) != d * nobs(obj)) stop("Clustering variable is not correct length.")  
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# nobs()

#-------------------------------------
# residuals
#-------------------------------------

residuals_CS.mlm <- function(obj) {
  res <- residuals(obj)
  as.vector(t(res))
}

#-------------------------------------
# model_matrix()
#-------------------------------------

model_matrix.mlm <- function(obj) {
  X <- model.matrix(obj)
  d <- ncol(residuals(obj))
  X_mat <- X %x% diag(1L, nrow = d)
  dimnames(X_mat) <- list(rep(dimnames(X)[[1]], each = d),
                          rep(dimnames(X)[[2]], times = d))
  X_mat
}

#----------------------------------------------
# get "working" variance-covariance matrix
#----------------------------------------------

targetVariance.mlm <- function(obj, cluster) {
  matrix_list(rep(1, nobs(obj) * ncol(residuals(obj))), cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.mlm <- function(obj, cluster) {
  weights <- weights(obj)
  if (is.null(weights)) {
    weights <- w_scale <- 1
  } else {
    w_scale <- mean(weights)
    weights <- weights / w_scale
  }
  W <- rep(weights, length.out = nobs(obj))
  W <- rep(W, each = ncol(residuals(obj)))
  W_list <- matrix_list(W, cluster, "both")
  attr(W_list, "w_scale") <- w_scale
  W_list
}

#----------------------------------------------
# get coefficient estimates
#----------------------------------------------

coef_CS.mlm <- function(obj) {
  as.vector(coef(obj))
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.mlm() is in sandwich package

v_scale.mlm <- function(obj) {
  nobs(obj)
}
