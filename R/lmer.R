#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an lmerMod object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from \code{\link[lme4:merMod-class]{merMod}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to 
#'   \code{getGroups(obj)}.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the
#'   estimated variance-covariance structure of the \code{lmerMod} object.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#'  
#' @examples 
#' library(lme4)
#' sleep_fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' vcovCR(sleep_fit, type = "CR2")
#' 
#' data(egsingle, package = "mlmRev")
#' math_model <- lmer(math ~ year * size + female + black + hispanic 
#'                    + (1 | schoolid) + (1 | childid), 
#'                    data = egsingle)
#' vcovCR(math_model, type = "CR2")
#' 
#' @export

vcovCR.lmerMod <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  
  if (!is.null(obj@call$weights)) stop("Models with prior weights are not currently supported.")
  
  if (missing(cluster)) cluster <- get_outer_group(obj)
  if(!is_nested_lmerMod(obj, cluster)) stop("Non-nested random effects detected. clubSandwich methods are not available for such models.")
  
  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

#-------------------------------------
# check nesting of random effects 
#-------------------------------------

get_outer_group <- function(obj) {
  group_n <- lme4::getME(obj, "l_i")
  group_facs <- lme4::getME(obj, "flist")
  group_facs[[which.min(group_n)]]
}

check_nested <- function(inner_grp, outer_grp) {
  n_outer <- tapply(outer_grp, inner_grp, function(x) length(unique(x)))
  all(n_outer == 1)
}

is_nested_lmerMod <- function(obj, cluster = get_outer_group(obj)) {
  group_facs <- lme4::getME(obj, "flist")
  nested <- vapply(group_facs, check_nested, outer_grp = cluster, FUN.VALUE = TRUE)
  all(nested)
}


# nobs()
# model_matrix()

#-------------------------------------
# coef_CS()
#-------------------------------------

coef_CS.lmerMod <- function(obj)
  lme4::getME(obj, "fixef")

#-------------------------------------
# residuals_CS()
#-------------------------------------

residuals_CS.lmerMod <- function(obj) {
  y <- lme4::getME(obj, "y")
  X <- lme4::getME(obj, "X")
  beta <- lme4::getME(obj, "fixef")
  y - as.numeric(X %*% beta)
}
  
  
#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.lmerMod <- function(obj, cluster = get_outer_group(obj)) {
  
  Z_mat <- lme4::getME(obj, "Z")
  Lambdat <- lme4::getME(obj, "Lambdat")
  Zlam_list <- matrix_list(Matrix::tcrossprod(Z_mat, Lambdat), fac = cluster, dim = "row")
  target_list <- lapply(Zlam_list, function(z) Matrix::tcrossprod(z) + Matrix::Diagonal(n = NROW(z)))
  
  lapply(target_list, as.matrix)
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.lmerMod <- function(obj, cluster = get_outer_group(obj)) {
  V_list <- targetVariance(obj, cluster)
  lapply(V_list, function(v) chol2inv(chol(v)))
  
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.lmerMod <- function(x, ...) {
  sigma_sq <- lme4::getME(x, "sigma")^2
  as.matrix(vcov(x) * v_scale(x)) / sigma_sq
}

v_scale.lmerMod <- function(obj) {
  min(lme4::getME(obj, "l_i"))
}
