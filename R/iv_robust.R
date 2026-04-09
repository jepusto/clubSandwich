#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an
#' \code{estimatr::iv_robust} object.
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates from an
#' \code{\link[estimatr]{iv_robust}} object.
#'
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{iv_robust} objects unless the model
#'   was fitted with a \code{clusters} argument.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be an identity matrix.
#' @param inverse_var Not used for \code{iv_robust} objects.
#' @inheritParams vcovCR
#'
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists
#'   of a matrix of the estimated variance of and covariances between the
#'   regression coefficient estimates.
#'
#' @seealso \code{\link{vcovCR}}
#'
#' @export

vcovCR.iv_robust <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ...) {

  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (missing(type)) stop("You must specify a type of sandwich estimator.")
  if (inverse_var != FALSE) stop("The inverse_var option is not available for iv_robust models.")

  obj$model.frame <- model.frame(obj)

  vcov_CR(obj, cluster = cluster, type = type,
          target = target, inverse_var = inverse_var, form = form)
}


#----------------------------------------------
# model.frame
#----------------------------------------------

#' @export

model.frame.iv_robust <- function(formula, ...) {

  mf <- formula$model.frame
  if (!is.null(mf)) return(mf)

  fit_env <- environment(formula$terms)

  cl <- formula$call

  mf_args <- match(c("formula", "data", "weights", "subset", "clusters"), names(cl), 0L)

  mf_cl <- cl[c(1L, mf_args)]
  mf_cl$formula <- formula$terms
  mf_cl[[1L]] <- quote(stats::model.frame)
  eval(mf_cl, envir = fit_env)
}


#----------------------------------------------
# get X matrix (projected)
#----------------------------------------------

#' @export

model_matrix.iv_robust <- function(obj) {

  mf <- model.frame(obj)

  # Regressor matrix X
  X <- model.matrix(obj$terms_regressors, data = mf)

  # Instrument formula from the RHS of the | in the IV formula
  inst_formula <- as.formula(paste("~", deparse(obj$formula[[3]][[3]])))
  Z <- model.matrix(inst_formula, data = mf)

  # Compute projected X: P_Z * X where P_Z = Z(Z'WZ)^{-1}Z'W
  w <- weights(obj)
  if (is.null(w)) {
    X_proj <- Z %*% solve(crossprod(Z)) %*% crossprod(Z, X)
  } else {
    X_proj <- Z %*% solve(crossprod(Z, w * Z)) %*% crossprod(Z, w * X)
    pos_wts <- w > 0
    if (!all(pos_wts)) {
      X_proj <- X_proj[pos_wts, , drop = FALSE]
    }
  }

  X_proj
}


#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.iv_robust <- function(x, ...) {

  N <- nobs(x)

  XZ <- model_matrix(x)

  if (x$weighted) {
    w <- weights(x)
    if (!is.null(w)) {
      pos_wts <- w > 0
      if (!all(pos_wts)) w <- w[pos_wts]
    }
    XtWX <- crossprod(XZ, w * XZ)
  } else {
    XtWX <- crossprod(XZ)
  }

  N * solve(XtWX)
}

#----------------------------------------------
# get residuals
#----------------------------------------------

#' @export

residuals_CS.iv_robust <- function(obj) {

  mf <- model.frame(obj)
  r <- mf[[obj$outcome]] - obj$fitted.values

  w <- weights(obj)
  if (!is.null(w)) {
    pos_wts <- w > 0
    if (!all(pos_wts)) r <- r[pos_wts]
  }

  r
}

# coef_CS() - use default
# targetVariance() - use default
# weightMatrix() - use default
# v_scale() - use default
