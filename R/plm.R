# coef()

#-----------------------------------------------
# Model matrix
#-----------------------------------------------
model_matrix.plm <- function(obj) {
  if (obj$args$model=="random") {
    model.matrix(Formula::as.Formula(formula(obj)), model.frame(obj))  
  } else {
    model.matrix(obj)
  }
}

#-------------------------------------
# unadjusted residuals
#-------------------------------------

residuals_CR.plm <- function(obj) {
  if (obj$args$model=="random") {
    y <- plm::pmodel.response(formula(obj), model.frame(obj), model = "pooling")
    Xb <- as.numeric(model_matrix(obj) %*% coef(obj))
    y - Xb
  } else {
    residuals(obj)
  }
}

#-------------------------------------
# nobs
#-------------------------------------

nobs.plm <- function(object, ...) {
  NROW(object$model)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.plm <- function(obj) {
  if (obj$args$model=="random") {
    if (obj$args$effect=="twoway") stop("Target variance is not block diagonal.")
    ind <- switch(obj$args$effect,
            individual = attr(model.frame(obj), "index")[[1]],
            time = attr(model.frame(obj), "index")[[1]])
    Z <- model.matrix(~ ind - 1)
    V <- tcrossprod(sqrt(obj$ercomp$sigma2[[3]]) * Z)
    diag(V) <- diag(V) + obj$ercomp$sigma2[[2]]
  } else {
    V <- rep(1, nobs(obj))
  }
  V
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.plm <- function(obj) {
  if (obj$args$model=="random") {
    if (obj$args$effect=="twoway") stop("Weight matrix is not block diagonal.")
    ind <- switch(obj$args$effect,
                  individual = attr(model.frame(obj), "index")[[1]],
                  time = attr(model.frame(obj), "index")[[1]])
    Z <- model.matrix(~ ind - 1)
    n_j <- colSums(Z)
    sigma_sq <- obj$ercomp$sigma2[[2]]
    tau_sq <- obj$ercomp$sigma2[[3]]
    theta_j <- tau_sq / ((n_j * tau_sq + sigma_sq) * sigma_sq)
    W <- diag(1 / sigma_sq, nrow = nobs(obj)) -  Z %*% (theta_j * t(Z))
  } else {
    W <- rep(1, nobs(obj))
  }
  W
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

vcovCR.plm <- function(obj, cluster, type, target, inverse_var) {
  
  if (obj$args$model=="random" & obj$args$effect=="twoway") stop("Variance matrix is not block diagonal.")
  
  if (missing(cluster)) {
    if (obj$args$effect=="twoway") stop("You must specify a clustering variable.")
    index <- attr(model.frame(obj),"index")
    cluster <- switch(obj$args$effect,
                      individual = index[[1]],
                      time = index[[2]])
  } else if (is.character(cluster)) {
    if (cluster %in% c("individual","time")) {
      index <- attr(model.frame(obj),"index")
      cluster <- switch(cluster,
                        individual = index[[1]],
                        time = index[[2]])
    }
  }

  if (missing(target)) target <- NULL
  if (missing(inverse_var) ) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}
