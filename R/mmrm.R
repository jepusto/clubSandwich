# Default methods work for mmrm objects:
#   model_matrix.default() — model.matrix(obj) works; weights are all positive
#   residuals_CS.default() — residuals(obj) works; weights are all positive
#   coef_CS.default() — coef(obj) works directly

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' @export

vcovCR.mmrm <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) {
    ff <- mmrm::component(obj, "full_frame")
    subject_var <- obj$formula_parts$subject_var
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
  vc <- VarCorr(obj)

  obs_by_subject <- split(seq_along(cluster), cluster)

  if (is.null(group_var)) {
    # Non-grouped: single covariance matrix for all subjects
    lapply(obs_by_subject, function(rows) {
      idx <- match(as.character(visits[rows]), visit_levels)
      unname(vc[idx, idx, drop = FALSE])
    })
  } else {
    # Grouped: pick the group-specific covariance matrix per subject
    groups <- ff[[group_var]]
    lapply(obs_by_subject, function(rows) {
      subj_group <- as.character(groups[rows[1]])
      idx <- match(as.character(visits[rows]), visit_levels)
      unname(vc[[subj_group]][idx, idx, drop = FALSE])
    })
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
