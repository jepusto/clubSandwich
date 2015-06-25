# coef()
# residuals_CR()
# vcov()
# model_matrix()

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.rma.uni <- function(obj) {
  obj$vi + obj$tau2
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.rma.uni <- function(obj) {
  if (obj$weighted) {
    if (is.null(obj$weights)) {
      wi <- 1 / (obj$vi + obj$tau2)  
    } else {
      wi <- obj$weights
    }
  } else {
    wi <- rep(1, obj$k)
  }
  wi
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.rma.uni <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (missing(target)) {
    target <- NULL
    inverse_var <- is.null(obj$weights) & obj$weighted
  } else {
    if (missing(inverse_var)) inverse_var <- FALSE
  }
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}
