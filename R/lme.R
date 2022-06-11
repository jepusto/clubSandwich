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
#' @examples 
#' 
#' if (requireNamespace("nlme", quietly = TRUE)) {
#' 
#'   library(nlme)
#'   rat_weight <- lme(weight ~ Time * Diet, data=BodyWeight, ~ Time | Rat) 
#'   vcovCR(rat_weight, type = "CR2")
#' 
#' }
#' 
#' if (requireNamespace("nlme", quietly = TRUE) & requireNamespace("mlmRev", quietly = TRUE)) {
#' 
#'   data(egsingle, package = "mlmRev")
#'   subset_ids <- levels(egsingle$schoolid)[1:10]
#'   math_model <- lme(math ~ year * size + female + black + hispanic, 
#'                     random = list(~ year | schoolid, ~ 1 | childid), 
#'                     data = egsingle, subset = schoolid %in% subset_ids)
#'   vcovCR(math_model, type = "CR2")
#'   
#' }
#' 
#' @export

vcovCR.lme <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
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

#' @export

residuals_CS.lme <- function(obj) 
  residuals(obj, level = 0)

#-------------------------------------
# coef_CS()
#-------------------------------------

#' @export

coef_CS.lme <- function(obj)
  nlme::fixef(obj)

#-------------------------------------
# model_matrix()
#-------------------------------------

#' @export

model_matrix.lme <- function(obj) {
  dat <- droplevels(getData(obj))
  model.matrix(formula(obj), data = dat)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

get_cor_grouping <- function(obj, levels = NULL) {
  if (!is.null(obj$groups)) {
    struct <- obj$modelStruct$corStruct
    if (is.null(struct)) struct <- obj
    mod_formula <- nlme::getGroupsFormula(struct)
    grps <- stats::model.frame(mod_formula, data = nlme::getData(obj))
    grps <- apply(grps, 1, paste, collapse = "/")
    if (is.null(levels)) levels <- unique(grps)
    grps <- factor(grps, levels = levels)
  } else if (!is.null(obj$modelStruct$corStruct)) {
    grps <- factor(rep("A",obj$dims$N))
  } else {
    grps <- factor(1:obj$dims$N)
  }
  grps
}

# Construct list of block-diagonal lowest-level var-cov matrices

get_sort_order <- function(obj) {
  groups <- obj$groups
  if (is.data.frame(groups)) {
    order(do.call(order, groups))
  } else if (!is.null(groups)) {
    order(order(groups))
  } else {
    1:obj$dims$N
  }
}

build_var_cor_mats <- function(obj) {
  
  if (is.null(obj$modelStruct$corStruct)) {
    
    # if there is no correlation structure,
    # then build block-diagonals with first available grouping variable
    
    if (is.null(obj$groups)) {
      
      # if there are no groups then make diagonal matrix-lists
      
      if (is.null(obj$modelStruct$varStruct)) {
        V_list <- as.list(rep(1, obj$dims$N))
      } else {
        sd_vec <- 1 / as.numeric(nlme::varWeights(obj$modelStruct$varStruct))
        V_list <- as.list(sd_vec^2)
      }
      grps <- factor(1:obj$dims$N)
      attr(V_list, "groups") <- grps
      names(V_list) <- levels(grps)
      
    } else {
      
      # if there are groups then make block-diagonal matrix-lists
      
      if (is.null(obj$modelStruct$varStruct)) {
        grps <- obj$groups[[1]]
        V_list <- tapply(rep(1, length(grps)),  grps, diag)
      } else {
        sort_order <- get_sort_order(obj)
        sd_vec <- 1 / as.numeric(nlme::varWeights(obj$modelStruct$varStruct))[sort_order]
        V_list <- tapply(sd_vec^2, obj$groups[[1]], diag)
      }
      attr(V_list, "groups") <- obj$groups[[1]]
    }
    
  } else {
    
    # if there is a correlation structure,
    # build block-diagonals according to its grouping structure
    
    R_list <- nlme::corMatrix(obj$modelStruct$corStruct)
    
    if (!is.list(R_list)) R_list <- list(A = R_list)
    
    grps <- get_cor_grouping(obj)
    missing_grps <- setdiff(levels(grps), names(R_list))
    
    if (length(missing_grps) > 0) {
      R_full <- rep(list(matrix(1,1,1)), nlevels(grps))
      names(R_full) <- levels(grps)
      R_full[names(R_list)] <- R_list
      R_list <- R_full
    } else {
      R_list <- R_list[levels(grps)]
    }
    
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- R_list
    } else {
      sort_order <- get_sort_order(obj)
      sd_vec <- 1 / as.numeric(nlme::varWeights(obj$modelStruct$varStruct))[sort_order]
      sd_list <- split(sd_vec, grps)
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    }
    
    attr(V_list, "groups") <- grps
  }

  return(V_list)
}

build_RE_mats <- function(obj) {
  
  # Get random effects structure
  all_groups <- rev(obj$groups)
  
  if (length(all_groups) == 1) {
    
    D_mat <- as.matrix(obj$modelStruct$reStruct[[1]])
    Z_mat <- model.matrix(obj$modelStruct$reStruc, nlme::getData(obj))
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)
    
    attr(ZDZ_list, "groups") <- all_groups[[1]]
    
  } else {
    D_list <- lapply(obj$modelStruct$reStruct, as.matrix)
    Z_mat <- model.matrix(obj$modelStruct$reStruc, nlme::getData(obj))
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")
    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)
    
    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]],
                                  big_mats = ZDZ_lists[[i]],
                                  crosswalk = all_groups[c(i-1,i)])
    }
    
    ZDZ_list <- ZDZ_lists[[i]]
    
    attr(ZDZ_list, "groups") <- all_groups[[i]]
    
  }
  
  ZDZ_list
  
}

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

#' @export

targetVariance.lme <- function(obj, cluster = nlme::getGroups(obj, level = 1)) {
  
  if (inherits(obj, "nlme")) stop("not implemented for \"nlme\" objects")
  
  # lowest-level covariance structure
  V_list <- build_var_cor_mats(obj)
  
  # random effects covariance structure
  ZDZ_list <- build_RE_mats(obj)
  
  V_grps <- attr(V_list, "groups")
  
  # Check if lowest-level covariance structure is nested within RE structure
  ZDZ_grps <- attr(ZDZ_list, "groups")
  group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
  nested <- all(group_mapping == 1L)
  
  if (nested) {
    target_list <- add_bdiag(V_list, ZDZ_list, data.frame(V_grps, ZDZ_grps))
    target_grps <- attr(ZDZ_list, "groups")
  } else {
    V_mat <- unblock(V_list, block = V_grps)
    ZDZ_mat <- unblock(ZDZ_list, block = ZDZ_grps)
    target_list <- V_mat + ZDZ_mat
    target_grps <- factor(rep("A", nrow(target_list)))
  }
  
  # check if clustering level is higher than highest level of random effects
  
  tb_groups <- table(target_grps)
  tb_cluster <- table(cluster)
  
  if (length(tb_groups) < length(tb_cluster) | 
      any(as.vector(tb_groups) != rep(as.vector(tb_cluster), length.out = length(tb_groups))) | 
      any(names(tb_groups) != rep(names(tb_cluster), length.out = length(tb_groups)))) {
    
    # check that random effects are nested within clusters  
    tb_cross <- table(target_grps, cluster)
    nested <- apply(tb_cross, 1, function(x) sum(x > 0) == 1)
    if (!all(nested)) stop("Random effects are not nested within clustering variable.")
    
    # expand target_list to level of clustering
    crosswalk <- data.frame(target_grps, cluster)
    target_list <- add_bdiag(small_mats = target_list, 
                             big_mats = matrix_list(rep(0, length(cluster)), cluster, dim = "both"),
                             crosswalk = crosswalk)
  }
  
  return(target_list)
}


targetVariance_old.lme <- function(obj, cluster = nlme::getGroups(obj, level = 1)) {
  
  if (inherits(obj, "nlme")) stop("not implemented for \"nlme\" objects")
  
  all_groups <- rev(obj$groups)
  smallest_groups <- all_groups[[1]]
  largest_groups <- all_groups[[length(all_groups)]]
  
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
    target_list <- Map("+", ZDZ_list, V_list)
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
    target_list <- ZDZ_lists[[i]]
  }
  
  # check if clustering level is higher than highest level of random effects
  
  tb_groups <- table(largest_groups)
  tb_cluster <- table(cluster)
  if (length(tb_groups) < length(tb_cluster) | 
      any(as.vector(tb_groups) != rep(as.vector(tb_cluster), length.out = length(tb_groups))) | 
      any(names(tb_groups) != rep(names(tb_cluster), length.out = length(tb_groups)))) {
    
    # check that random effects are nested within clusters  
    tb_cross <- table(largest_groups, cluster)
    nested <- apply(tb_cross, 1, function(x) sum(x > 0) == 1)
    if (!all(nested)) stop("Random effects are not nested within clustering variable.")
    
    # expand target_list to level of clustering
    crosswalk <- data.frame(largest_groups, cluster)
    target_list <- add_bdiag(small_mats = target_list, 
                             big_mats = matrix_list(rep(0, length(cluster)), cluster, dim = "both"),
                             crosswalk = crosswalk)
  }
  
  return(target_list)
}


#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weightMatrix.lme <- function(obj, cluster = nlme::getGroups(obj, level = 1)) {
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

#' @export

v_scale.lme <- function(obj) {
  nlevels(nlme::getGroups(obj))
}
