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
#'   Clustering on individuals can also be obtained by specifying the name of
#'   the individual index (e.g., \code{cluster = "state"}) or \code{cluster =
#'   "individual"}; clustering on time periods can be obtained by specifying the
#'   name of the time index (e.g., \code{cluster = "year"}) or \code{cluster =
#'   "time"}; if a group index is specified, clustering on groups (in which
#'   individuals are nested) can be obtained by specifying the name of the group
#'   index or \code{cluster = "group"}. For random-effects models, the cluster
#'   will be taken equal to the included random effect identifier if not
#'   otherwise specified.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. By default, the target is taken to be an identity
#'   matrix for fixed effect models or the estimated compound-symmetric
#'   covariance matrix for random effects models.
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
#'               data = Produc, index = c("state","year","region"),
#'               effect = "individual", model = "within")
#' vcovCR(plm_FE, type="CR2")
#' vcovCR(plm_FE, type = "CR2", cluster = Produc$region) # clustering on region
#'
#' # random effects
#' plm_RE <- update(plm_FE, model = "random")
#' vcovCR(plm_RE, type = "CR2")
#' vcovCR(plm_RE, type = "CR2", cluster = Produc$region) # clustering on region
#'
#' # nested random effects
#' plm_nested <- update(plm_FE, effect = "nested", model = "random")
#' vcovCR(plm_nested, type = "CR2") # clustering on region
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
  
  if (inherits(dat, "pdata.frame") | is.numeric(index_names)) {
    indices <- plm::index(obj)
  } else {
    if (is.null(index_names)) index_names <- names(dat)[1:2]
    indices <- as.list(dat[index_names])
  }
  
  if (length(indices) == 3) indices <- indices[c(3,1:2)]
  do.call(order, args = indices)
}

findCluster.plm <- function(obj, cluster) {
  
  index <- attr(model.frame(obj),"index")
  effect <- obj$args$effect
  
  if (obj$args$model=="random" & effect=="twoways") stop("Variance matrix is not block diagonal. clubSandwich methods are not available for such models.")
  
  if (missing(cluster)) {
    
    # Infer missing clustering variable
    
    if (effect=="twoways") stop("You must specify a clustering variable when effect = 'twoways'.")
    
    cluster <- switch(obj$args$effect,
                      individual = index[[1]],
                      time = index[[2]],
                      nested = index[[3]])
    
  } else if ((length(cluster)==1) & is.character(cluster)) {
    
    # Translate clustering variable keywords
    
    allowed_clusters <- switch(effect, 
                               individual = c("individual","group"),
                               time = "time",
                               twoways = "none",
                               nested = "group")
    
    if (cluster %in% c("individual","time","group")) {
      
      # Check for nesting of random effects inside clusters
      if (obj$args$model == "random" & !cluster %in% allowed_clusters) {
          err_msg <- paste0("For a random effects model, cluster = '", cluster, "' is not allowed when effect = '", effect, "'.")
          stop(err_msg)
      }
      
      cluster <- switch(cluster,
                        individual = index[[1]],
                        time = index[[2]],
                        group = index[[3]])
      
    } else if (cluster %in% names(index)) {

      cluster_ID <- c("individual","time","group")[1:length(index)][cluster == names(index)]
      
      # Check for nesting of random effects inside clusters
      if (obj$args$model == "random" & !cluster_ID %in% allowed_clusters) {
        err_msg <- paste0("For a random effects model, cluster = '", cluster, "' is not allowed when effect = '", effect, "'.")
        stop(err_msg)
      }
      
      cluster <- index[[cluster]]
            
    } else {
      err_msg <- paste0("Clustering variable could not be inferred. Please check the argument cluster = '", cluster, "'.")
      stop(err_msg)
    }
    
  } else {
    
    # Sort user-specified clustering variable
    
    sort_order <- get_index_order(obj)
    cluster <- cluster[sort_order]
    
  }

  # Check for nesting of random effects inside clusters
  
  if (obj$args$model == "random") {
    RE_index <- switch(effect, 
                       individual = index[[1]],
                       time = index[[2]],
                       nested = index[[3]])
    
    if(!check_nested(RE_index, cluster)) stop("Random effects are not nested within clustering variable. clubSandwich methods are not available for such models.")
  }
  
  if (obj$args$model=="fd") {
    cluster <- cluster[index[[2]] != levels(index[[2]])[1]]
  }
  
  cluster
}

#-----------------------------------------------
# Model matrix
#-----------------------------------------------

#' @export

model_matrix.plm <- function(obj) {
  if (obj$args$model=="random") {
    model.matrix(Formula::as.Formula(formula(obj)), model.frame(obj))  
  } else {
    cstcovar.rm <- switch(obj$args$model,
                          within = "all",
                          fd = "covariates",
                          pooling = "none",
                          between = "none")
    model.matrix(obj, 
                 model = obj$args$model, 
                 effect = obj$args$effect,
                 cstcovar.rm = cstcovar.rm)
  }
}

#----------------------------------------------
# Augmented model matrix
#----------------------------------------------

#' @export

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

#' @export

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

block_mat <- function(nj, r) {
  Vj <- matrix(r, nj, nj)
  diag(Vj) <- 1 + r
  Vj
}

#' @export

targetVariance.plm <- function(obj, cluster) {
  
  if (obj$args$model=="random") {
    
    index <- attr(model.frame(obj),"index")
    r <- obj$ercomp$sigma2[[2]] / obj$ercomp$sigma2[[1]]
    
    if (obj$args$effect=="nested") {
      
      r_grp <- obj$ercomp$sigma2[[3]] / obj$ercomp$sigma2[[1]]
      grp_mat <- lapply(table(index[[3]]), function(x) matrix(r_grp, nrow = x, ncol = x))
      ind_mat <- lapply(table(index[[1]]), block_mat, r = r)
      RE_index <- index[[3]]
      target_mat <- add_bdiag(ind_mat, grp_mat, crosswalk = index[c(1,3)])
      
    } else {
    
      RE_index <- switch(obj$args$effect,
                         individual = index[[1]],
                         time = index[[2]])
      target_mat <- lapply(table(RE_index), block_mat, r = r)
      
    }
    
    nest_bdiag(target_mat, crosswalk = data.frame(RE_index, cluster))
    
  } else {
    
    matrix_list(rep(1, nobs(obj)), cluster, "both")
    
  }
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

inverse_block_mat <- function(nj, r) {
  theta <- r / (1 + r * nj)
  Wj <- matrix(-theta, nj, nj)
  diag(Wj) <- 1 - theta
  Wj
}

inverse_nested_block_mat <- function(nj, r1, r2) {
  V_inv <- lapply(nj, inverse_block_mat, r = r1)
  const <- 1 / sqrt(1 / r2 + sum(nj / (1 + r1 * nj)))
  vec <- const * rep(1 / (1 + nj * r1), nj)
  indices <- factor(names(vec), levels = names(nj))
  add_submatrices(indices, small_mat = V_inv, big_mat = -tcrossprod(vec))
}

#' @export

weightMatrix.plm <- function(obj, cluster) {
  
  if (obj$args$model=="random") {
    
    index <- attr(model.frame(obj),"index")
    r <- obj$ercomp$sigma2[[2]] / obj$ercomp$sigma2[[1]]
    
    if (obj$args$effect=="nested") {
      
      r_grp <- obj$ercomp$sigma2[[3]] / obj$ercomp$sigma2[[1]]
      njs <- tapply(index[[1]], index[[3]], function(x) table(droplevels(x)))
      RE_index <- index[[3]]
      target_mat <- lapply(njs, inverse_nested_block_mat, r1 = r, r2 = r_grp)
      
    } else {
      
      RE_index <- switch(obj$args$effect,
                         individual = index[[1]],
                         time = index[[2]])
      target_mat <- lapply(table(RE_index), inverse_block_mat, r = r)
      
    }
    
    nest_bdiag(target_mat, crosswalk = data.frame(RE_index, cluster))
    
  } else {
    
    matrix_list(rep(1, nobs(obj)), cluster, "both")
    
  }
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.plm <- function(x, ...) {
  v_scale(x) * vcov(x) / with(x, sum(residuals^2) / df.residual)
}

#' @export

v_scale.plm <- function(obj) {
  max(sapply(attr(obj$model, "index"), nlevels))
}
