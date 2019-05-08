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
#' @param ignore_FE Optional logical controlling whether fixed effects are
#'   ignored when calculating small-sample adjustments in models where fixed
#'   effects are estimated through absorption.
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
#' library(plm)
#' # fixed effects
#' data("Produc", package = "plm")
#' plm_FE <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#'               data = Produc, index = c("state","year"), 
#'               effect = "individual", model = "within")  
#' vcovCR(plm_FE, type="CR2")
#' 
#' # random effects
#' plm_RE <- update(plm_FE, model = "random")
#' vcovCR(plm_RE, type = "CR2")
#' 
#' # first differencing
#' data(Fatalities, package = "AER")
#' Fatalities <- within(Fatalities, {
#'   frate <- 10000 * fatal / pop
#'   drinkagec <- cut(drinkage, breaks = 18:22, include.lowest = TRUE, right = FALSE)
#'   drinkagec <- relevel(drinkagec, ref = 4)
#' })
#' 
#' plm_FD <- plm(frate ~ beertax + drinkagec + miles + unemp + log(income), 
#'               data = Fatalities, index = c("state", "year"), 
#'               model = "fd")
#' vcovHC(plm_FD, method="arellano", type = "sss", cluster = "group")
#' vcovCR(plm_FD, type = "CR1S")
#' vcovCR(plm_FD, type = "CR2")
#' 
#' 
#' @export

vcovCR.plm <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ignore_FE = FALSE, ...) {
  
  if (obj$args$model=="random" & obj$args$effect=="twoways") stop("Variance matrix is not block diagonal.")
  
  if (missing(cluster)) {
    cluster <- findCluster.plm(obj = obj)
  } else {
    cluster <- findCluster.plm(obj = obj, cluster = cluster)  
  } 
  
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  obj$na.action <- attr(obj$model, "na.action")
  
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, 
          form = form, ignore_FE = ignore_FE)
}

get_index_order <- function(obj) {
  envir <- environment(obj$formula)
  mf <- match.call(plm::plm, call = obj$call, envir = envir)
  dat <- eval(mf$data, envir)
  
  index_names <- eval(mf$index)
  
  if ("pdata.frame" %in% class(dat) | is.numeric(index_names)) {
    indices <- index(obj)
  } else {
    if (is.null(index_names)) index_names <- names(dat)[1:2]
    indices <- as.list(dat[index_names])
  }
  
  do.call(order, args = indices)
}

findCluster.plm <- function(obj, cluster) {
  index <- attr(model.frame(obj),"index")
  if (missing(cluster)) {
    if (obj$args$effect=="twoways") stop("You must specify a clustering variable.")
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
  
  if (obj$args$model=="fd") {
    cluster <- cluster[index[[2]] != levels(index[[2]])[1]]
  }
  
  cluster
}

#-----------------------------------------------
# Model matrix
#-----------------------------------------------

model_matrix.plm <- function(obj) {
  if (obj$args$model=="random") {
    model.matrix(Formula::as.Formula(formula(obj)), model.frame(obj))  
  } else {
    model.matrix(obj, model = obj$args$model, effect = obj$args$effect)
  }
}

#----------------------------------------------
# Augmented model matrix
#----------------------------------------------

augmented_model_matrix.plm <- function(obj, cluster, inverse_var, ignore_FE) {
  index <- attr(model.frame(obj),"index")
  individual <- droplevels(as.factor(index[[1]]))
  time <- droplevels(as.factor(index[[2]]))
  effect <- obj$args$effect
  
  if (ignore_FE) {
    S <- NULL 
  } else if (obj$args$model=="within") {
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
        S <- residuals(lm.fit(model.matrix(~ 0 + individual), model.matrix(~ 0 + time)[,-1]))
      } else if (inverse_var & identical(time, cluster)) {
        S <- residuals(lm.fit(model.matrix(~ 0 + time), model.matrix(~ 0 + individual)[,-1]))
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

residuals_CS.plm <- function(obj) {
  if (obj$args$model=="random") {
    y <- plm::pmodel.response(formula(obj), model.frame(obj), model = "pooling")
    nm <- names(y)
    y <- as.numeric(y)
    names(y) <- nm
    Xb <- as.numeric(model_matrix(obj) %*% coef(obj))
    y - Xb
  } else {
    residuals(obj)
  }
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.plm <- function(obj, cluster) {
  if (obj$args$model=="random") {
    block_mat <- function(nj) {
      sigma_sq <- obj$ercomp$sigma2[[1]]
      tau_sq <- obj$ercomp$sigma2[[2]]
      r <- tau_sq / sigma_sq
      Vj <- matrix(r, nj, nj)
      diag(Vj) <- 1 + r
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
    sigma_sq <- obj$ercomp$sigma2[[1]]
    tau_sq <- obj$ercomp$sigma2[[2]]
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

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

bread.plm <- function(x, ...) {
  # if (x$args$model=="random") {
  #   v_scale(x) * vcov(x) / x$ercomp$sigma2$idios
  # } else {
  #   v_scale(x) * vcov(x) / with(x, sum(residuals^2) / df.residual) 
  # }
  v_scale(x) * vcov(x) / with(x, sum(residuals^2) / df.residual)
}

v_scale.plm <- function(obj) {
  max(sapply(attr(obj$model, "index"), nlevels))
}
