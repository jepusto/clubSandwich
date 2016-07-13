#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a plm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a \code{\link[plm]{plm}} 
#' object.
#' 
#' @param cluster Optional character string, expression, or vector indicating 
#'   which observations belong to the same cluster. For fixed-effect models that
#'   include individual effects or time effects (but not both), the cluster will
#'   be taken equal to the included fixed effects if not otherwise specified. 
#'   Clustering on individuals can also be obtained by taking \code{cluster = 
#'   "individual"} and clustering on time periods can be obtained with 
#'   \code{cluster = "time"}. For random-effects models, the cluster will be 
#'   taken equal to the included random effect identifier if not otherwise 
#'   specified.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. By default, the target is taken to be an identity
#'   matrix for fixed effect models or the estimated compound-symmetric covariance matrix for random effects models. 
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

vcovCR.plm <- function(obj, cluster, type, target, inverse_var) {
  
  if (obj$args$model=="random" & obj$args$effect=="twoways") stop("Variance matrix is not block diagonal.")
  
  index <- attr(model.frame(obj),"index")
  
  if (missing(cluster)) {
    if (obj$args$effect=="twoways") stop("You must specify a clustering variable.")
    index <- attr(model.frame(obj),"index")
    cluster <- switch(obj$args$effect,
                      individual = index[[1]],
                      time = index[[2]])
  } else if ((length(cluster)==1) & is.character(cluster)) {
    if (cluster %in% c("individual","time")) {
      cluster <- switch(cluster,
                        individual = index[[1]],
                        time = index[[2]])
    } 
  } else {
    sort_order <- get_index_order(obj)
    cluster <- cluster[sort_order]
  }
  
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  obj$na.action <- attr(obj$model, "na.action")
  
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

get_index_order <- function(obj) {
  envir <- environment(obj$formula)
  mf <- match.call(plm::plm, call = obj$call, envir = envir)
  index_names <- names(attr(model.frame(obj), "index"))
  index <- eval(mf$data, envir)[,index_names]
  order(index[,1],index[,2])
}

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

#----------------------------------------------
# Augmented model matrix
#----------------------------------------------

augmented_model_matrix.plm <- function(obj, cluster, inverse_var) {
  index <- attr(model.frame(obj),"index")
  individual <- droplevels(as.factor(index[[1]]))
  time <- droplevels(as.factor(index[[2]]))
  effect <- obj$args$effect
  
  if (obj$args$model=="within") {
    if (effect=="individual") {
      if (inverse_var & identical(individual, cluster)) {
        S <- NULL
      } else {
        S <- model.matrix(~ 0 + individual)
      }
    } else if (effect=="time") {
      if (inverse_var & identical(time, cluster)) {
        S <- NULL 
      } else {
        S <- model.matrix(~ 0 + time)
      }
    } else if (effect=="twoways") {
      if (inverse_var & identical(individual, cluster)) {
        S <- model.matrix(~ 0 + time)
      } else if (inverse_var & identical(time, cluster)) {
        S <- model.matrix(~ 0 + individual)
      } else {
        S <- model.matrix(~ 0 + individual + time)
      }
    } 
  } else {
    S <- NULL
  }
  
  return(S)
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

targetVariance.plm <- function(obj, cluster) {
  if (obj$args$model=="random") {
    block_mat <- function(nj) {
      Vj <- matrix(obj$ercomp$sigma2$id, nj, nj)
      diag(Vj) <- obj$ercomp$sigma2$idios + obj$ercomp$sigma2$id
      Vj
    }
    lapply(table(cluster), block_mat)
  } else {
    matrix_list(rep(1, nobs(obj)), cluster, "both")
  }
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.plm <- function(obj, cluster) {
  if (obj$args$model=="random") {
    sigma_sq <- obj$ercomp$sigma2$idios
    tau_sq <- obj$ercomp$sigma2$id
    block_mat <- function(nj) {
      theta <- tau_sq / ((nj * tau_sq + sigma_sq))
      Wj <- matrix(-theta, nj, nj)
      diag(Wj) <- 1 - theta
      Wj
    }
    lapply(table(cluster), block_mat)
  } else {
    matrix_list(rep(1, nobs(obj)), cluster, "both")
  }
}
