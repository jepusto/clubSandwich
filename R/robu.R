
#-----------------------------------------------
# coefficients
#-----------------------------------------------

coef.robu <- function(object, ...) {
  beta <- as.vector(object$b.r)
  labs <- object$reg_table$labels
  names(beta) <- levels(labs)[labs]
  beta
}

#-----------------------------------------------
# residuals
#-----------------------------------------------

residuals.robu <- function(object, ...) {
  ord <- order(order(object$study_orig_id))
  object$data.full$e.r[ord]
}


#-----------------------------------------------
# Model matrix
#-----------------------------------------------

model.matrix.robu <- function(object, ...) {
  ord <- order(order(object$study_orig_id))
  object$Xreg[ord,]
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  1 / obj$data.full$r.weights[ord]
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.robu <- function(obj) {
  ord <- order(order(obj$study_orig_id))
  if (obj$user_weighting) { 
    obj$data.full$userweights[ord]
  } else{
    obj$data.full$r.weights[ord]
  }
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.robu <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- obj$study_orig_id
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}
