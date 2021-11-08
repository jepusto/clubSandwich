#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a glm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{\link{glm}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{glm} objects.
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
#' data(dietox, package = "geepack")
#' dietox$Cu <- as.factor(dietox$Cu)
#' weight_fit <- glm(Weight ~ Cu * poly(Time, 3), data=dietox, family = "quasipoisson")
#' V_CR <- vcovCR(weight_fit, cluster = dietox$Pig, type = "CR2")
#' coef_test(weight_fit, vcov = V_CR, test = "Satterthwaite")
#' 
#' @export

vcovCR.glm <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
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

model_matrix.glm <- function(obj) {
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

residuals_CS.glm <- function(obj) {
  residuals(obj, type = "response")
}

#-----------------------------------------------
# Get (model-based) working variance matrix 
#-----------------------------------------------

#' @export

targetVariance.glm <- function(obj, cluster) {
  mu <- fitted.values(obj)
  var_fun <- obj$family$variance
  v <- var_fun(mu)
  w <- weights(obj, type = "prior")
  matrix_list(v / w, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weightMatrix.glm <- function(obj, cluster) {
  mu <- fitted.values(obj)
  var_fun <- obj$family$variance
  v <- var_fun(mu)
  w <- weights(obj, type = "prior")
  matrix_list(w / v, cluster, "both")
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

# bread.glm() is in sandwich package

#' @export

v_scale.glm <- function(obj) {
  if (substr(obj$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
    dispersion <- 1
  } else {
    wres <- as.vector(residuals(obj, "working")) * weights(obj, "working")
    dispersion <- sum(wres^2)/sum(weights(obj, "working"))
  } 
  as.vector(sum(summary(obj)$df[1:2])) * dispersion
}
