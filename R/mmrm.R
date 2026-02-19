#-------------------------------------
# model_matrix()
#-------------------------------------

#' @export

model_matrix.mmrm = function(obj) {
  model.matrix(obj)
}

#-------------------------------------
# residuals_CS()
#-------------------------------------

#' @export

residuals_CS.mmrm = function(obj) {
  residuals(obj)
}

#-------------------------------------
# coef_CS()
#-------------------------------------

#' @export

coef_CS.mmrm = function(obj) {
  coef(obj)
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' @export

vcovCR.mmrm = function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) {
    ff = mmrm::component(obj, "full_frame")
    subject_var = obj$formula_parts$subject_var
    cluster = droplevels(as.factor(ff[[subject_var]]))
  }
  if (missing(target)) target = NULL
  if (missing(inverse_var)) inverse_var = is.null(target)
  vcov_CR(obj, cluster = cluster, type = type,
          target = target, inverse_var = inverse_var, form = form)
}

#-------------------------------------
# targetVariance()
#  Returns a list of per-subject working covariance matrices (Psi_j)
#  For non-grouped models, obj$cov is the full visit x visit matrix where each subject's block is the submatrix for their observed visits
#  For grouped models, obj$cov is a named list (one matrix per group)
#-------------------------------------

# < Code function here > 

#-------------------------------------
# weightMatrix()
#  Returns the inverse (precision) of each subjects covariance block
#-------------------------------------


# < Code function here > 
