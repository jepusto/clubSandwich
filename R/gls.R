#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a gls object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from a \code{\link[nlme]{gls}} object.
#' 
#' @param cluster Optional expression or vector indicating which observations 
#'   belong to the same cluster. If not specified, will be set to 
#'   \code{getGroups(obj)}.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If not specified, the target is taken to be the
#'   estimated variance-covariance structure of the \code{gls} object.
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
#' library(nlme)
#' data(Ovary, package = "nlme")
#' Ovary$time_int <- 1:nrow(Ovary)
#' lm_AR1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary, 
#'               correlation = corAR1(form = ~ time_int | Mare))
#' vcovCR(lm_AR1, type = "CR2")
#'     
#' @export

vcovCR.gls <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  if (missing(cluster)) cluster <- nlme::getGroups(obj)
  if (missing(target)) target <- NULL
  if (missing(inverse_var) ) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# residuals_CS()
# coef()
# nobs()

#-------------------------------------
# model_matrix()
#-------------------------------------

getData <- function (object) {
  if ("data" %in% names(object)) {
    data <- object$data
  } else {
    dat_name <- object$call$data
    envir_names <- sys.frames()
    ind <- sapply(envir_names, function(e) exists(as.character(dat_name), envir = e))
    e <- envir_names[[min(which(ind))]]
    data <- eval(dat_name, envir = e)
  }
  if (is.null(data)) return(data)
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    data <- if (inherits(naAct, "omit")) {
      data[-naAct, ]
      
    } else if (inherits(naAct, "exclude")) {
      data
    } else eval(object$call$na.action)(data)
  }
  subset <- object$call$subset
  if (!is.null(subset)) {
    subset <- eval(asOneSidedFormula(subset)[[2]], data)
    data <- data[subset, ]
  }
  data
}

#' @export

model_matrix.gls <- function(obj) {
  dat <- getData(obj)
  model.matrix(formula(obj), data = dat)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

#' @export

targetVariance.gls <- function(obj, cluster = nlme::getGroups(obj)) {
  
  groups <- nlme::getGroups(obj)
  if (is.null(groups)) groups <- cluster
  
  if (is.null(obj$modelStruct$corStruct)) {
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- matrix_list(rep(1, length(cluster)), cluster, "both")
    } else {
      wts <- nlme::varWeights(obj$modelStruct$varStruct)
      V_list <- matrix_list(1 / wts^2, cluster, "both")
    } 
  } else {
    R_list <- nlme::corMatrix(obj$modelStruct$corStruct)
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- R_list
    } else {
      sd_vec <- 1 / nlme::varWeights(obj$modelStruct$varStruct)[order(order(groups))]
      sd_list <- split(sd_vec, groups)
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    } 
  } 
  
  # check if clustering level is higher than highest level of random effects
  
  tb_groups <- table(groups)
  tb_cluster <- table(cluster)
  if (length(tb_groups) < length(tb_cluster) | 
      any(as.vector(tb_groups) != rep(as.vector(tb_cluster), length.out = length(tb_groups))) | 
      any(names(tb_groups) != rep(names(tb_cluster), length.out = length(tb_groups)))) {
    
    # check that random effects are nested within clusters  
    tb_cross <- table(groups, cluster)
    nested <- apply(tb_cross, 1, function(x) sum(x > 0) == 1)
    if (!all(nested)) stop("Random effects are not nested within clustering variable.")
    
    # expand target_list to level of clustering
    crosswalk <- data.frame(groups, cluster)
    V_list <- add_bdiag(small_mats = V_list, 
                             big_mats = matrix_list(rep(0, length(cluster)), cluster, dim = "both"),
                             crosswalk = crosswalk)
  }
  
  V_list
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

#' @export

weightMatrix.gls <- function(obj, cluster = nlme::getGroups(obj)) {
  V_list <- targetVariance(obj, cluster)
  lapply(V_list, function(v) chol2inv(chol(v)))
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.gls <- function(x, ...) {
  vcov(x) * nobs(x) / x$sigma^2
}

# v_scale() is default