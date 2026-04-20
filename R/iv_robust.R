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

  if (inverse_var != FALSE) stop("The inverse_var option is not available for iv_robust models.")

  if (isTRUE(obj$fes) && !requireNamespace("fixest", quietly = TRUE)) message("For improved performance in models with fixed effects, install the package {fixest}.")

  obj$model.frame <- model.frame(obj)

  if (missing(cluster)) {
    cluster <- findCluster.iv_robust(obj)
    if (is.null(cluster)) stop("You must specify a clustering variable or `obj` must include a clustering variable.")
  }

  if (missing(type)) {
    type <- switch(obj$se_type, CR0 = "CR0", CR2 = "CR2", stata = "CR1S", "No valid SE type")
    if (type == "No valid SE type") stop("You must specify a `type` of sandwich estimator to calculate or `obj` must include an `se_type` of 'CR0', 'CR2', or 'stata'.")
  }

  vcov_CR(obj, cluster = cluster, type = type,
          target = target, inverse_var = inverse_var, form = form)
}


#' Pulls the clustering variable from an iv_robust object, if it has one.
#'
#' @param obj an iv_robust object
#'
#' @return The data within the clustering variable
#' @keywords internal
#' @noRd

findCluster.iv_robust <- function(obj) {

  if (!obj$clustered) return(NULL)

  model.frame(obj)[["(clusters)"]]
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
  mf <- eval(mf_cl, envir = fit_env)

  # If no fixed effects, we're done
  if (!isTRUE(formula$fes)) return(mf)

  # Build a second model.frame for the fixed-effects variables
  fe_args <- match(c("fixed_effects", "data", "subset"), names(cl), 0L)
  fe_cl <- cl[c(1L, fe_args)]
  names(fe_cl)[[2]] <- "formula"
  fe_cl[[2]] <- reformulate(all.vars(fe_cl[[2]]))
  fe_cl[[1L]] <- quote(stats::model.frame)
  mf_fe <- eval(fe_cl, envir = fit_env)

  # Reconcile any differences in omitted rows between the two model frames
  mf_omit <- na.action(mf)
  if (!is.null(names(mf_omit))) mf_omit <- names(mf_omit)
  fe_omit <- na.action(mf_fe)
  if (!is.null(names(fe_omit))) fe_omit <- names(fe_omit)

  if (identical(mf_omit, fe_omit)) {
    mf_combined <- cbind(mf, mf_fe)
    mf_combined_omit <- mf_omit
  } else {
    mf_combined <- cbind(
      mf[!(rownames(mf) %in% fe_omit), , drop = FALSE],
      mf_fe[!(rownames(mf_fe) %in% mf_omit), , drop = FALSE]
    )
    mf_combined_omit <- sort(c(mf_omit, fe_omit))
    i_unique <- !duplicated(mf_combined_omit)
    mf_combined_omit <- mf_combined_omit[i_unique]
    class(mf_combined_omit) <- "omit"
  }

  attr(mf_combined, "terms") <- attr(mf, "terms")
  attr(mf_combined, "na.action") <- mf_combined_omit

  return(mf_combined)
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

  w <- weights(obj)

  # With fixed effects, residualize X and Z on the FE design matrix
  # (Frisch-Waugh-Lovell); the intercept is absorbed by the FEs.
  if (isTRUE(obj$fes)) {

    intercept_X <- colnames(X) == "(Intercept)"
    if (any(intercept_X)) X <- X[, !intercept_X, drop = FALSE]
    intercept_Z <- colnames(Z) == "(Intercept)"
    if (any(intercept_Z)) Z <- Z[, !intercept_Z, drop = FALSE]

    fe_formula <- as.formula(obj$call$fixed_effects)

    if (requireNamespace("fixest", quietly = TRUE)) {
      frame <- mf[attr(terms(fe_formula), "term.labels")]
      X <- fixest::demean(X = X, f = frame, weights = w)
      Z <- fixest::demean(X = Z, f = frame, weights = w)
    } else {
      fe_formula <- update(fe_formula, ~ . - 1)
      varnames <- all.vars(fe_formula)
      for (v in varnames) mf[[v]] <- as.factor(mf[[v]])
      F_mat <- model.matrix(fe_formula, data = mf)
      if (is.null(w)) {
        X <- stats::lm.fit(F_mat, X)$residuals
        Z <- stats::lm.fit(F_mat, Z)$residuals
      } else {
        X <- stats::lm.wfit(F_mat, X, w)$residuals
        Z <- stats::lm.wfit(F_mat, Z, w)$residuals
      }
    }
  }

  # Compute projected X: X_proj = Z (Z'WZ)^{-1} Z'W X
  if (is.null(w)) {
    X_proj <- Z %*% solve(crossprod(Z), crossprod(Z, X))
  } else {
    X_proj <- Z %*% solve(crossprod(Z, w * Z), crossprod(Z, w * X))
    pos_wts <- w > 0
    if (!all(pos_wts)) {
      X_proj <- X_proj[pos_wts, , drop = FALSE]
    }
  }

  X_proj
}


#----------------------------------------------
# augmented model matrix (fixed effects)
#----------------------------------------------

#' @export

augmented_model_matrix.iv_robust <- function(obj, cluster, inverse_var, ignore_FE) {

  if (!isTRUE(obj$fes)) return(NULL)

  fe_formula <- as.formula(obj$call$fixed_effects)
  fe_formula <- update(fe_formula, ~ . - 1)

  mf <- model.frame(obj)

  varnames <- all.vars(fe_formula)
  for (v in varnames) mf[[v]] <- as.factor(mf[[v]])

  model.matrix(fe_formula, data = mf)
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

#----------------------------------------------
# na.action
#----------------------------------------------

#' @export

na.action.iv_robust <- function(object, ...) {
  na.action(model.frame(object))
}

# coef_CS() - use default
# targetVariance() - use default
# weightMatrix() - use default
# v_scale() - use default
