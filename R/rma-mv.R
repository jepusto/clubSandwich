
#----------------------------------------------------------------------
# utility function for computing block-diagonal covariance matrices
#----------------------------------------------------------------------

#' Impute a block-diagonal covariance matrix
#'
#' @description \loadmathjax{} \code{impute_covariance_matrix} calculates a
#'   block-diagonal covariance matrix, given the marginal variances, the block
#'   structure, and an assumed correlation structure. Can be used to create
#'   compound-symmetric structures, AR(1) auto-correlated structures, or
#'   combinations thereof.
#'
#' @param vi Vector of variances
#' @param cluster Vector indicating which effects belong to the same cluster.
#'   Effects with the same value of `cluster` will be treated as correlated.
#' @param r Vector or numeric value of assumed constant correlation(s) between
#'   effect size estimates from each study.
#' @param ti Vector of time-points describing temporal spacing of effects, for
#'   use with auto-regressive correlation structures.
#' @param ar1 Vector or numeric value of assumed AR(1) auto-correlation(s)
#'   between effect size estimates from each study. If specified, then \code{ti}
#'   argument must be specified.
#' @param smooth_vi Logical indicating whether to smooth the marginal variances
#'   by taking the average \code{vi} within each cluster. Defaults to
#'   \code{FALSE}.
#' @param subgroup Vector of category labels describing sub-groups of effects.
#'   If non-null, effects that share the same category label and the same
#'   cluster will be treated as correlated, but effects with different category
#'   labels will be treated as uncorrelated, even if they come from the same
#'   cluster.
#' @param return_list Optional logical indicating whether to return a list of
#'   matrices (with one entry per block) or the full variance-covariance matrix.
#'
#'
#'
#' @return If \code{cluster} is appropriately sorted, then a list of matrices,
#'   with one entry per cluster, will be returned by default. If \code{cluster}
#'   is out of order, then the full variance-covariate matrix will be returned
#'   by default. The output structure can be controlled with the optional
#'   \code{return_list} argument.
#'
#' @details A block-diagonal variance-covariance matrix (possibly represented as a list
#'   of matrices) with a specified structure. The structure depends on whether
#'   the \code{r} argument, \code{ar1} argument, or both arguments are
#'   specified. Let \mjeqn{v_{ij}}{v-ij} denote the specified variance for
#'   effect \mjeqn{i}{i} in cluster \mjeqn{j}{j} and \mjeqn{C_{hij}}{C-hij} be
#'   the covariance between effects \mjeqn{h}{h} and \mjeqn{i}{i} in cluster
#'   \mjeqn{j}{j}.  
#'   \itemize{ 
#'     \item{If only \code{r} is specified,}{ each block of
#'     the variance-covariance matrix will have a constant (compound symmetric)
#'     correlation, so that \mjdeqn{C_{hij} = r_j \sqrt{v_{hj} v_{ij},}}{C-hij = r-j *
#'     sqrt(v-hj v-ij),} where \mjeqn{r_j}{r-j} is the specified correlation for cluster \mjeqn{j}{j}. 
#'     If only a single value is given in \code{r}, then it will be used for every cluster.} 
#'     \item{If only \code{ar1} is specified,}{ each block of the
#'     variance-covariance matrix will have an AR(1) auto-correlation structure,
#'     so that \mjdeqn{C_{hij} = \phi_j^{|t_{hj} - t_{ij}|} \sqrt{v_{hj}
#'     v_{ij},}}{C-hij = (ar1-j)^|t-hj - t-ij| * sqrt(v-hj v-ij),} where
#'     where \mjeqn{\phi_j}{ar1-j} is the specified auto-correlation for cluster \mjeqn{j}{j} and 
#'     \mjeqn{t_{hj}}{t-hj} and \mjeqn{t_{ij}}{t-ij} are specified time-points
#'     corresponding to effects \mjeqn{h}{h} and \mjeqn{i}{i} in cluster
#'     \mjeqn{j}{j}. If only a single value is given in \code{ar1}, then it will be used 
#'     for every cluster.}
#'     \item{If both \code{r} and \code{ar1} are specified,}{ each block of the
#'     variance-covariance matrix will have combination of compound symmetric and 
#'     an AR(1) auto-correlation structures,
#'     so that \mjdeqn{C_{hij} = \left[r_j + (1 - r_j)\phi_j^{|t_{hj} - t_{ij}|}\right] \sqrt{v_{hj}
#'     v_{ij},}}{C-hij = [r-j + (1 - r-j)(ar1-j)^|t-hj - t-ij|] * sqrt(v-hj v-ij),} where
#'     where where \mjeqn{r_j}{r-j} is the specified constant correlation for cluster \mjeqn{j}{j}, 
#'     \mjeqn{\phi_j}{ar1-j} is the specified auto-correlation for cluster \mjeqn{j}{j} and 
#'     \mjeqn{t_{hj}}{t-hj} and \mjeqn{t_{ij}}{t-ij} are specified time-points
#'     corresponding to effects \mjeqn{h}{h} and \mjeqn{i}{i} in cluster
#'     \mjeqn{j}{j}. If only single values are given in \code{r} or \code{ar1}, they will be used 
#'     for every cluster.}
#'   }
#'   If \code{smooth_vi = TRUE},
#'   
#' @export
#'
#' @examples
#' library(metafor)
#' 
#' # Constant correlation
#' data(SATcoaching)
#' V_list <- impute_covariance_matrix(vi = SATcoaching$V, cluster = SATcoaching$study, r = 0.66)
#' MVFE <- rma.mv(d ~ 0 + test, V = V_list, data = SATcoaching)
#' conf_int(MVFE, vcov = "CR2", cluster = SATcoaching$study)
#' 


impute_covariance_matrix <- function(vi, cluster, r, ti, ar1, 
                                     smooth_vi = FALSE, 
                                     subgroup = NULL, 
                                     return_list = identical(as.factor(cluster), sort(as.factor(cluster)))) {
  
  cluster <- droplevels(as.factor(cluster))
  
  vi_list <- split(vi, cluster)
  
  if (smooth_vi) vi_list <- lapply(vi_list, function(x) rep(mean(x, na.rm = TRUE), length(x)))

  if (missing(r) & missing(ar1)) stop("You must specify a value for r or for ar1.")
  
  if (!missing(r)) {
    r_list <- rep_len(r, length(vi_list))
    if (missing(ar1)) {
      vcov_list <- Map(function(V, rho) (rho + diag(1 - rho, nrow = length(V))) * tcrossprod(sqrt(V)), 
                       V = vi_list, 
                       rho = r_list)
    }
  } 
  
  if (!missing(ar1)) {
    if (missing(ti)) stop("If you specify a value for ar1, you must provide a vector for ti.")
    
    ti_list <- split(ti, cluster)
    ar_list <- rep_len(ar1, length(vi_list))
    
    if (missing(r)) {
      vcov_list <- Map(function(V, time, phi) (phi^as.matrix(stats::dist(time))) * tcrossprod(sqrt(V)), 
                       V = vi_list, 
                       time = ti_list, 
                       phi = ar_list)
    } else {
      vcov_list <- Map(function(V, rho, time, phi) (rho + (1 - rho) * phi^as.matrix(stats::dist(time))) * tcrossprod(sqrt(V)), 
                       V = vi_list, 
                       rho = r_list, 
                       time = ti_list, 
                       phi = ar_list)
    }
    
    vcov_list <- lapply(vcov_list, function(x) {
      attr(x, "dimnames") <- NULL
      x
    })
  } 
  
  if (!is.null(subgroup)) {
    si_list <- split(subgroup, cluster)
    subgroup_list <- lapply(si_list, function(x) sapply(x, function(y) y == x))
    vcov_list <- Map(function(V, S) V * S, V = vcov_list, S = subgroup_list)
  }
  
  if (return_list) {
    return(vcov_list)
  } else {
    vcov_mat <- metafor::bldiag(vcov_list)
    cluster_index <- order(order(cluster))
    return(vcov_mat[cluster_index, cluster_index])
  }
}

#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a robu object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a 
#' \code{\link[metafor]{rma.mv}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to the factor in
#'   the random-effects structure with the fewest distinct levels. Caveat
#'   emptor: the function does not check that the random effects are nested.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the 
#'   estimated variance-covariance structure of the \code{rma.mv} object.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export
#' 
#' @examples
#' library(metafor)
#' data(hierdat, package = "robumeta")
#' 
#' mfor_fit <- rma.mv(effectsize ~ binge + followup + sreport + age, 
#'                  V = var, random = list(~ 1 | esid, ~ 1 | studyid),
#'                  data = hierdat)
#' mfor_fit
#' 
#' mfor_CR2 <- vcovCR(mfor_fit, type = "CR2")
#' mfor_CR2
#' coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
#' 
#' Wald_test(mfor_fit, constraints = constrain_zero(c(2,4)), vcov = mfor_CR2)
#' Wald_test(mfor_fit, constraints = constrain_zero(2:5), vcov = mfor_CR2)

vcovCR.rma.mv <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  
  if (obj$withR) stop("vcovCR.rma.mv() does not work with fixed correlation matrices in the R argument.")
  
  if (missing(cluster)) {
    cluster <- findCluster.rma.mv(obj)
  } else {
    # check that random effects are nested within clustering variable
    mod_struct <- parse_structure(obj)
    
    if (length(cluster) != NROW(mod_struct$cluster_dat)) {
      cluster <- cluster[obj$not.na]
    } 
    
    nested <- test_nested(cluster, fac = mod_struct$cluster_dat)
    if (!all(nested)) stop("Random effects are not nested within clustering")
  }
  
  if (missing(target)) {
    target <- NULL
    inverse_var <- is.null(obj$W)
  } else {
    if (missing(inverse_var)) inverse_var <- FALSE
  }
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# coef()
# residuals_CS()
# vcov()
# model_matrix

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.rma.mv <- function(obj, cluster) {
  matrix_list(obj$M, cluster, "both")
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.rma.mv <- function(obj, cluster) {
  if (is.null(obj$W)) {
    V_list <- targetVariance(obj, cluster)
    lapply(V_list, function(v) chol2inv(chol(v)))
  } else{
    matrix_list(obj$W, cluster, "both")
  }
}

#-----------------------------------------------
# Get outer-most clustering variable
#-----------------------------------------------

get_structure <- function(obj) {
  data.frame(G = obj$withG, H = obj$withH, R = obj$withR, S = obj$withS)
}

test_nested <- function(cluster, fac) {
  
  if (is.list(fac)) {
    res <- sapply(fac, test_nested, cluster = cluster)
    return(res)
  } 
  
  groupings <- tapply(cluster, fac, function(x) length(unique(x)))
  all(groupings==1L)  
}

parse_structure <- function(obj) {
  
  level_dat <- vector(mode = "integer")
  cluster_dat <- data.frame(row.names = 1:obj$k)
  
  if (obj$withG) {
    level_dat[["G"]] <- obj$g.nlevels[[2]]
    cluster_dat$G <- obj$mf.g$outer
  }
  
  if (obj$withH) {
    level_dat[["H"]] <- obj$h.nlevels[[2]]
    cluster_dat$H <- obj$mf.h$outer
  }
  
  if (obj$withS) {
    s_levels <- obj$s.nlevels
    names(s_levels) <- obj$s.names
    level_dat <- c(level_dat, s_levels)
    
    mf_all <- do.call(cbind, obj$mf.r)
    mf_s <- mf_all[obj$s.names]
    cluster_dat <- cbind(cluster_dat, mf_s)
  }
  
  list(level_dat = level_dat, cluster_dat = cluster_dat)
}

findCluster.rma.mv <- function(obj) {
  
  if (obj$withR) stop("vcovCR.rma.mv() does not work with fixed correlation matrices in the R argument.")
  
  # parse model structure
  mod_struct <- parse_structure(obj) 
  
  if (length(mod_struct$level_dat) == 0L) stop("No clustering variable specified.")
  
  # determine cluster with smallest number of levels
  
  highest_cluster <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  cluster <- mod_struct$cluster_dat[[highest_cluster]]
  
  # check that random effects are nested within clustering variable
  nested <- test_nested(cluster, fac = mod_struct$cluster_dat)
  if (!all(nested)) stop("Random effects are not nested within clustering")
  
  # clean up
  if (!is.factor(cluster)) cluster <- as.factor(cluster)
  droplevels(cluster)
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

bread.rma.mv <- function(x, ...) {
  if (is.null(x$W)) {
    B <- vcov(x) * nobs(x)
  } else{
    X_mat <- model_matrix(x)
    XWX <- t(X_mat) %*% x$W %*% X_mat
    B <- chol2inv(chol(XWX)) * nobs(x)
    rownames(B) <- colnames(B) <- colnames(X_mat)
  }
  B
}

v_scale.rma.mv <- function(obj) {
  nobs(obj)
}
