# Default methods that work for fixest objects:
#   coef_CS.default() — coef(obj) returns covariate coefficients (no FE)
#   residuals_CS.default() — residuals(obj) returns full residuals (including FE)
#   targetVariance.default() — identity target is appropriate for OLS
#   weightMatrix.default() — weights(obj) works correctly
#   v_scale.default() — nobs(obj) works correctly

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a fixest object.
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates from a
#' \code{\link[fixest]{feols}} object.
#'
#' @param cluster Optional expression or vector indicating which observations
#'   belong to the same cluster. If not specified, will be set to the fixed
#'   effect variable when the model has exactly one fixed effect; for models
#'   with multiple fixed effects, the clustering variable must be specified
#'   explicitly.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be an identity matrix.
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
#' if (requireNamespace("fixest", quietly = TRUE)) withAutoprint({
#'
#'   library(fixest)
#'   data("ChickWeight", package = "datasets")
#'   ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)
#'
#'   # Single fixed effect: cluster auto-detected from the FE
#'   fit_one <- feols(weight ~ Time | Chick, data = ChickWeight)
#'   vcovCR(fit_one, type = "CR2")
#'   coef_test(fit_one, vcov = "CR2", test = "Satterthwaite")
#'
#'   # Two-way fixed effects: cluster must be specified explicitly
#'   fit_two <- feols(weight ~ Time | Chick + Diet, data = ChickWeight)
#'   vcovCR(fit_two, cluster = ChickWeight$Chick, type = "CR2")
#'
#' })
#'
#' @export

vcovCR.fixest <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) {
    cluster <- findCluster.fixest(obj)
    if (is.null(cluster)) stop("You must specify a clustering variable.")
  }
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type,
          target = target, inverse_var = inverse_var, form = form)
}


#' Pulls clustering variable from fixest objects.
#'
#' @param obj a fixest object
#'
#' @return A factor vector of cluster membership, or NULL
#' @keywords internal
#' @noRd

findCluster.fixest <- function(obj) {
  if (length(obj$fixef_vars) != 1L) return(NULL)

  fe_var <- obj$fixef_vars

  fe_mf <- tryCatch(model.matrix(obj, type = "fixef"),
                    error = function(e) NULL)
  if (is.null(fe_mf) || is.null(fe_mf[[fe_var]])) return(NULL)

  fe_mf[[fe_var]]
}


#-------------------------------------
# model_matrix()
#-------------------------------------

#' @export

model_matrix.fixest <- function(obj) {

  # If feols() was called with demeaned = TRUE, use the stored demeaned matrix
  if (length(obj$fixef_vars) > 0L && !is.null(obj$X_demeaned)) {
    X <- obj$X_demeaned
  } else if (length(obj$fixef_vars) > 0L) {
    X <- model.matrix(obj)
    fe_df <- data.frame(lapply(obj$fixef_id, factor), check.names = FALSE)
    wts <- weights(obj)
    # Slope FEs (Destination[Year], Destination[[Year]]) need slope.vars and
    # slope.flag passed through, otherwise fixest::demean residualizes on the
    # intercept FE only and X is wrong.
    has_slopes <- !is.null(obj$slope_flag) && any(obj$slope_flag != 0)
    if (has_slopes && !is.null(wts)) {
      X <- fixest::demean(X = X, f = fe_df,
                          slope.vars = obj$slope_variables_reordered,
                          slope.flag = obj$slope_flag,
                          weights = wts)
    } else if (has_slopes) {
      X <- fixest::demean(X = X, f = fe_df,
                          slope.vars = obj$slope_variables_reordered,
                          slope.flag = obj$slope_flag)
    } else if (!is.null(wts)) {
      X <- fixest::demean(X = X, f = fe_df, weights = wts)
    } else {
      X <- fixest::demean(X = X, f = fe_df)
    }
  } else {
    X <- model.matrix(obj)
  }

  # Handle zero weights
  w <- weights(obj)
  if (!is.null(w)) {
    pos_wts <- w > 0
    if (!all(pos_wts)) {
      X <- X[pos_wts, , drop = FALSE]
    }
  }

  X
}

#-------------------------------------
# augmented_model_matrix()
#-------------------------------------

#' @export

augmented_model_matrix.fixest <- function(obj, cluster, inverse_var, ignore_FE) {
  if (ignore_FE || length(obj$fixef_vars) == 0L) return(NULL)

  # slope_flag entries: 0 = intercept-only FE, 1 = intercept + slope (D[Y]),
  # -1 = slope-only (D[[Y]]). For each FE we emit intercept dummies if the
  # flag is >= 0 and slope columns (dummy * slope_var) if the flag is non-zero.
  flags <- if (is.null(obj$slope_flag)) {
    rep(0L, length(obj$fixef_vars))
  } else {
    obj$slope_flag
  }
  # slope_variables_reordered is indexed positionally over FE terms with slopes
  # (those where flags != 0), in the same order as fixef_vars.
  slope_vars <- obj$slope_variables_reordered
  slope_idx <- 0L

  blocks <- list()
  for (i in seq_along(obj$fixef_vars)) {
    fe_name <- obj$fixef_vars[i]
    f <- factor(obj$fixef_id[[i]])
    dum <- model.matrix(~ 0 + f)

    if (flags[i] >= 0) {
      colnames(dum) <- paste0(fe_name, levels(f))
      blocks[[length(blocks) + 1L]] <- dum
    }
    if (flags[i] != 0) {
      slope_idx <- slope_idx + 1L
      sv <- slope_vars[[slope_idx]]
      slope_block <- dum * sv
      colnames(slope_block) <- paste0(fe_name, ":", levels(f))
      blocks[[length(blocks) + 1L]] <- slope_block
    }
  }

  do.call(cbind, blocks)
}

#---------------------------------------
# bread and scaling constant
#---------------------------------------

# bread.fixest() is in the fixest package
# v_scale.default() uses nobs(), which works for fixest
