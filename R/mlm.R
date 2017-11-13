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
  rownames(X_mat) <- rep(dimnames(X)[[1]], each = d)
  colnames(X_mat) <- paste(rep(colnames(residuals(obj)), ncol(X)), 
                           rep(colnames(X), each = d), sep = ":")
  i <- unlist(lapply(1:d, function(x) seq(x, ncol(X_mat), d)))
  X_mat[,i]
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
  cf <- coef(obj)
  res <- as.vector(cf)
  names(res) <- paste(rep(colnames(cf), each = nrow(cf)), 
                      rep(rownames(cf), ncol(cf)), sep = ":")
  res
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export
#' 

bread.mlm <- function (x, ...) {
  if (!is.null(x$na.action)) class(x$na.action) <- "omit"
  X_mat <- model.matrix(x)
  w <- weights(x)
  XWX <- if (!is.null(w)) crossprod(X_mat, w * X_mat) else crossprod(X_mat)
  B <- chol2inv(chol(XWX))
  rval <- diag(ncol(residuals(x))) %x% (B * nobs(x))
  
  col_names <- paste(rep(colnames(residuals(x)), each = ncol(X_mat)), 
                     rep(colnames(X_mat), ncol(residuals(x))), sep = ":")
  colnames(rval) <- rownames(rval) <- col_names
  return(rval)
}


v_scale.mlm <- function(obj) {
  nobs(obj)
}
