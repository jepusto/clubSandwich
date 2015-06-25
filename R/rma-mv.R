# coef()
# residuals_CR()
# vcov()
# model_matrix

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.rma.mv <- function(obj) {
  obj$M
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.rma.mv <- function(obj) {
  if (is.null(obj$W)) { 
    chol2inv(chol(obj$M))
  } else{
    obj$W
  }
}

#-----------------------------------------------
# Get outer-most clustering variable
#-----------------------------------------------

findCluster.rma.mv <- function(obj) {
  if (obj$withS) {
    r <- which.min(obj$s.nlevels)
    cluster <- obj$mf.r[[r]][[obj$s.names[r]]]
  } else if (obj$withG) {
    cluster <- obj$mf.r[[1]][[obj$g.names[2]]]
  } else {
    stop("No clustering variable specified.")
  }
  droplevels(as.factor(cluster))
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.rma.mv <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- findCluster.rma.mv(obj)
  if (missing(target)) {
    target <- NULL
    inverse_var <- is.null(obj$W)
  } else {
    if (missing(inverse_var)) inverse_var <- FALSE
  }
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}
