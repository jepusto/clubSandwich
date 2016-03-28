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
#' @export

vcovCR.gls <- function(obj, cluster, type, target, inverse_var) {
  if (missing(cluster)) cluster <- nlme::getGroups(obj)
  if (missing(target)) target <- NULL
  if (missing(inverse_var) ) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target, inverse_var = inverse_var)
}

# residuals_CR()
# coef()
# nobs()

#-------------------------------------
# model_matrix()
#-------------------------------------

getData <- function (object) {
  mCall <- object$call
  dat_name <- if ("data" %in% names(object)) object$data else mCall$data
  envir_names <- sys.frames()
  ind <- sapply(envir_names, function(e) exists(as.character(dat_name), envir = e))
  e <- envir_names[[min(which(ind))]]
  data <- eval(dat_name, envir = e)
  if (is.null(data)) return(data)
  naAct <- object[["na.action"]]
  if (!is.null(naAct)) {
    data <- if (inherits(naAct, "omit")) {
      data[-naAct, ]
      
    } else if (inherits(naAct, "exclude")) {
      data
    } else eval(mCall$na.action)(data)
  }
  subset <- mCall$subset
  if (!is.null(subset)) {
    subset <- eval(asOneSidedFormula(subset)[[2]], data)
    data <- data[subset, ]
  }
  data
}

model_matrix.gls <- function(obj) {
  dat <- getData(obj)
  model.matrix(formula(obj), data = dat)
}

#-------------------------------------
# Get (model-based) working variance matrix 
#-------------------------------------

targetVariance.gls <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  V <- matrix(0, N, N)
  for (i in levels(groups)) {
    V[groups == i, groups == i] <- nlme::getVarCov(obj, individual = i)
  }
  V
}

#-------------------------------------
# Get weighting matrix
#-------------------------------------

weightMatrix.gls <- function(obj) {
  groups <- nlme::getGroups(obj)
  N <- nobs(obj)
  W <- matrix(0, N, N)
  for (i in levels(groups)) {
    W[groups == i, groups == i] <- chol2inv(chol(nlme::getVarCov(obj, individual = i)))
  }
  W
}
