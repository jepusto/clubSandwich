# Default methods work for mmrm objects:
#   model_matrix.default() — model.matrix(obj) works; weights are all positive
#   residuals_CS.default() — residuals(obj) works; weights are all positive
#   coef_CS.default() — coef(obj) works directly

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an mmrm object.
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates from an
#' \code{\link[mmrm]{mmrm}} object.
#'
#' @param cluster Optional expression or vector indicating which observations
#'   belong to the same cluster. If not specified, will be set to the subject
#'   variable from the \code{mmrm} model formula.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If not specified, the target is taken to be the
#'   estimated variance-covariance structure of the \code{mmrm} object, with
#'   weights applied if present.
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
#' if (requireNamespace("mmrm", quietly = TRUE)) withAutoprint({
#'
#'   library(mmrm)
#'   data(fev_data, package = "mmrm")
#'
#'   # Fit an mmrm model with unstructured covariance
#'   mmrm_fit <- mmrm(
#'     FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
#'     data = fev_data
#'   )
#'
#'   # CR2 cluster-robust variance estimator (cluster auto-detected)
#'   vcovCR(mmrm_fit, type = "CR2")
#'
#'   # Coefficient tests with Satterthwaite degrees of freedom
#'   coef_test(mmrm_fit, vcov = "CR2", test = "Satterthwaite")
#'
#'   # Fit a weighted mmrm model
#'   fev_data$wt <- 1 + rpois(nrow(fev_data), lambda = 3)
#'   mmrm_wt <- mmrm(
#'     FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
#'     data = fev_data,
#'     weights = fev_data$wt
#'   )
#'
#'   # CR2 works with weighted models
#'   vcovCR(mmrm_wt, type = "CR2")
#'
#' })
#'
#' @export

vcovCR.mmrm <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) {
    ff <- mmrm::component(obj, "full_frame")
    subject_var <- mmrm::component(obj, "subject_var")
    cluster <- droplevels(as.factor(ff[[subject_var]]))
  }
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type,
          target = target, inverse_var = inverse_var, form = form)
}

#-------------------------------------
# targetVariance()
#-------------------------------------

#' @export

targetVariance.mmrm <- function(obj, cluster) {
  ff <- mmrm::component(obj, "full_frame")
  visit_var <- obj$formula_parts$visit_var
  group_var <- obj$formula_parts$group_var
  visit_levels <- levels(ff[[visit_var]])
  visits <- ff[[visit_var]]
  vc <- mmrm::VarCorr(obj)

  obs_by_subject <- split(seq_along(cluster), cluster)

  if (is.null(group_var)) {
    # Non-grouped: single covariance matrix for all subjects
    V_list <- lapply(obs_by_subject, function(rows) {
      idx <- match(as.character(visits[rows]), visit_levels)
      unname(vc[idx, idx, drop = FALSE])
    })
  } else {
    # Grouped: pick the group-specific covariance matrix per subject
    groups <- ff[[group_var]]
    V_list <- lapply(obs_by_subject, function(rows) {
      subj_group <- as.character(groups[rows[1]])
      idx <- match(as.character(visits[rows]), visit_levels)
      unname(vc[[subj_group]][idx, idx, drop = FALSE])
    })
  }

  # Multiply by weights if non-unit
  wts <- ff[["(weights)"]]

  if (all(wts == 1)) {
    return(V_list)
  } else {
    wt_inv_sqrt_j <- split(1 / sqrt(wts), cluster)
    V_weighted <- Map(function(iw, V) tcrossprod(iw) * V, iw = wt_inv_sqrt_j, V = V_list)
    return(V_weighted)
  }
}

#-------------------------------------
# weightMatrix()
#-------------------------------------

#' @export

weightMatrix.mmrm <- function(obj, cluster) {
  V_list <- targetVariance(obj, cluster)
  lapply(V_list, function(v) chol2inv(chol(v)))
}

#---------------------------------------
# bread and scaling constant
#---------------------------------------

#' @export

bread.mmrm <- function(x, ...) {
  vcov(x) * v_scale(x)
}

#' @export

v_scale.mmrm <- function(obj) {
  mmrm::component(obj, "n_subjects")
}

