#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a geeglm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{\link{geeglm}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{geeglm} objects.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be the estimated variance function.
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
#' if (requireNamespace("geepack", quietly = TRUE)) {
#' 
#'   data(dietox, package = "geepack")
#'   dietox$Cu <- as.factor(dietox$Cu)
#'   mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))
#'   gee1 <- geeglm(mf, data=dietox, id=Pig, family=poisson("identity"), corstr="ar1")
#'   V_CR <- vcovCR(gee1, cluster = dietox$Pig, type = "CR2")
#'   coef_test(gee1, vcov = V_CR, test = "Satterthwaite")
#'   
#' }
#' 
#' @export

vcovCR.geeglm <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# coef()
# nobs()

#-----------------------------------------------
# Model matrix
#-----------------------------------------------

#' @export

model_matrix.geeglm <- function(obj) {
  X <- model.matrix(obj)
  eta <- obj$linear.predictors
  dmu_deta <- obj$family$mu.eta
  d <- dmu_deta(eta)
  d * X
}

#-------------------------------------
# residuals
#-------------------------------------

#' @export

residuals_CS.geeglm <- function(obj) {
  residuals(obj, type = "response")
}

#-----------------------------------------------
# Get (model-based) working variance matrix 
#-----------------------------------------------

#' @export

targetVariance.geeglm <- function(obj, cluster) {
  mu <- fitted.values(obj)
  var_fun <- obj$family$variance
  v <- as.numeric(var_fun(mu))
  a <- tapply(v, obj$id, sqrt)
  aa <- lapply(a, Matrix::tcrossprod)
  if (obj$corstr %in% c("independence", "exchangeable", "ar1", "unstructured", "userdefined") == F) {
    stop("Working correlation matrix must be a matrix with the following correlation structures: independence, exchangeable, ar1, unstructured, or userdefined")
  } 
  else if 
     (obj$corstr == "ar1") {
      ar1_cor <- function(n, alpha) {
        exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                          (1:n - 1))
        alpha^exponent
      }
      r_ar1 <- lapply(obj$geese$clusz, ar1_cor, alpha = obj$geese$alpha)
      v <- mapply("*", aa, r_ar1)
      v <- unlist(v) # To solve
      w <- weights(obj, type = "prior")
      matrix_list(v / w, cluster, "both")
     }
  else {
    paste("0") # Try "ar1" for now, will add others
  }
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weightMatrix.geeglm <- function(obj, cluster) {
  mu <- fitted.values(obj)
  var_fun <- obj$family$variance
  v <- as.numeric(var_fun(mu))
  w <- weights(obj, type = "prior")
  matrix_list(w / v, cluster, "both")
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.glm() is in sandwich package

#' @export

v_scale.geeglm <- function(obj) {
  if (substr(obj$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
    dispersion <- 1
  } else {
    wres <- as.vector(residuals(obj, "working")) * weights(obj, "working")
    dispersion <- sum(wres^2)/sum(weights(obj, "working"))
  } 
  as.vector(sum(summary(obj)$df[1:2])) * dispersion
}
