#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an lme object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a \code{\link[nlme]{lme}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to 
#'   \code{getGroups(obj)}.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the
#'   estimated variance-covariance structure of the \code{lme} object.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

vcovCR.lme <- function(obj, cluster, type, target, inverse_var, form = "sandwich") {
  # if (length(obj$groups) > 1) stop("vcovCR.lme does not work for models with multiple levels of random effects.")
  if (missing(cluster)) cluster <- nlme::getGroups(obj, level = 1)
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# nobs()

#-------------------------------------
# residuals_CS()
#-------------------------------------

residuals_CS.lme <- function(obj) 
  residuals(obj, level = 0)

#-------------------------------------
# coef_CS()
#-------------------------------------

coef_CS.lme <- function(obj)
  nlme::fixef(obj)

#-------------------------------------
# model_matrix()
#-------------------------------------

model_matrix.lme <- function(obj) {
  dat <- droplevels(getData(obj))
  model.matrix(formula(obj), data = dat)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

targetVariance.lme <- function(obj, cluster) {
  
  if (any("nlme" == class(obj))) stop("not implemented for \"nlme\" objects")
  
  all_groups <- rev(obj$groups)
  # if (length(all_groups) > 1) stop("not implemented for multiple levels of nesting")
  smallest_groups <- all_groups[[1]]
  
  # Get level-1 variance-covariance structure as V_list
  
  if (is.null(obj$modelStruct$corStruct)) {
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- matrix_list(rep(1, length(smallest_groups)), smallest_groups, "both")
    } else {
      wts <- nlme::varWeights(obj$modelStruct$varStruct)[order(do.call(order, all_groups))]
      V_list <- matrix_list(1 / wts^2, smallest_groups, "both")
    } 
  } else {
    R_list <- as.list(rep(1, nlevels(smallest_groups)))
    names(R_list) <- levels(smallest_groups)
    R_sublist <- nlme::corMatrix(obj$modelStruct$corStruct)
    R_list[names(R_sublist)] <- R_sublist
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- R_list
    } else {
      sd_vec <- 1 / nlme::varWeights(obj$modelStruct$varStruct)[order(do.call(order, all_groups))]
      sd_list <- split(sd_vec, smallest_groups)
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    } 
  } 
  
  # Get random effects structure
  
  if (length(all_groups) == 1) {
    D_mat <- as.matrix(obj$modelStruct$reStruct[[1]])
    Z_mat <- model.matrix(obj$modelStruct$reStruc, getData(obj))
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)
    return(Map("+", ZDZ_list, V_list)      )
  } else {
    D_list <- lapply(obj$modelStruct$reStruct, as.matrix)
    Z_mat <- model.matrix(obj$modelStruct$reStruc, getData(obj))
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")
    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)
    ZDZ_lists[[1]] <- Map("+", ZDZ_lists[[1]], V_list)  
    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]], 
                                  big_mats = ZDZ_lists[[i]], 
                                  crosswalk = all_groups[c(i-1,i)])
    }
    return(ZDZ_lists[[i]])
  }
}


#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.lme <- function(obj, cluster) {
  V_list <- targetVariance(obj, cluster)
  lapply(V_list, function(v) chol2inv(chol(v)))
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.lme <- function(x, ...) {
  vcov(x) * v_scale(x) / x$sigma^2
}

v_scale.lme <- function(obj) {
  nlevels(nlme::getGroups(obj))
}
