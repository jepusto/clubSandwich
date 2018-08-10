
#----------------------------------------------------------------------
# utility function for computing block-diagonal covariance matrices
#----------------------------------------------------------------------

#' Impute a block-diagonal covariance matrix
#' 
#' \code{impute_covariance_matrix} calculates a block-diagonal covariance 
#' matrix, given the marginal variances, the block structure, and an assumed 
#' correlation.
#' 
#' @param vi Vector of variances
#' @param cluster Vector indicating which effects belong to the same cluster. 
#'   Effects with the same value of `cluster` will be treated as correlated.
#' @param r Vector or numeric value of assume correlation(s) between effect size
#'   estimates from each study.
#' @param return_list Optional logical indicating whether to return a list of matrices
#'   (with one entry per block) or the full variance-covariance matrix.
#'   
#' @return If \code{cluster} is appropriately sorted, then a list of matrices, 
#'   with one entry per cluster, will be returned by default. If \code{cluster}
#'   is out of order, then the full variance-covariate matrix will be returned
#'   by default. The output structure can be controlled with the optional
#'   \code{return_list} argument.
#'   
#' @export
#' 
#' @examples
#' library(metafor)
#' data(SATcoaching)
#' V_list <- impute_covariance_matrix(vi = SATcoaching$V, cluster = SATcoaching$study, r = 0.66)
#' MVFE <- rma.mv(d ~ 0 + test, V = V_list, data = SATcoaching)
#' coef_test(MVFE, vcov = "CR2", cluster = SATcoaching$study)
#' 


impute_covariance_matrix <- function(vi, cluster, r, return_list = identical(as.factor(cluster), sort(as.factor(cluster)))) {
  
  cluster <- droplevels(as.factor(cluster))
  
  vi_list <- split(vi, cluster)
  r_list <- rep_len(r, length(vi_list))
  vcov_list <- Map(function(V, rho) (rho + diag(1 - rho, nrow = length(V))) * tcrossprod(sqrt(V)), V = vi_list, rho = r_list)
  
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
#' Wald_test(mfor_fit, constraints = c(2,4), vcov = mfor_CR2)
#' Wald_test(mfor_fit, constraints = 2:5, vcov = mfor_CR2)

vcovCR.rma.mv <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) cluster <- findCluster.rma.mv(obj)
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
