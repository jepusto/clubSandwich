#----------------------------------------------
# user-facing vcovCR function
#----------------------------------------------

#' Cluster-robust variance-covariance matrix
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @param obj Fitted model for which to calcualte the variance-covariance matrix
#' @param cluster Expression or vector indicating which observations belong to 
#'   the same cluster. For some classes, the cluster will be detected 
#'   automatically if not specified.
#' @param type Character string specifying which small-sample adjustment should 
#'   be used.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If a vector, the target matrix is assumed to be 
#'   diagonal. If not specified, \code{vcovCR} will attempt to infer a value.
#' @param inverse_var Optional logical indicating whether the weights used in 
#'   fitting the model are inverse-variance. If not specified, \code{vcovCR} 
#'   will attempt to infer a value.
#' @param form Controls the form of the returned matrix. The default 
#'   \code{"sandwich"} will return the sandwich variance-covariance matrix. 
#'   Alternately, setting \code{form = "meat"} will return only the meat of the
#'   sandwich and setting \code{form = B}, where \code{B} is a matrix of
#'   appropriate dimension, will return the sandwich variance-covariance matrix
#'   calculated using \code{B} as the bread.
#'   
#' @description This is a generic function, with specific methods defined for 
#'   \code{\link[stats]{lm}}, \code{\link[plm]{plm}}, \code{\link[nlme]{gls}}, 
#'   \code{\link[nlme]{lme}}, \code{\link[robumeta]{robu}}, 
#'   \code{\link[metafor]{rma.uni}}, and \code{\link[metafor]{rma.mv}} objects.
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates. The matrix has several attributes: 
#'   \describe{ \item{type}{indicates which small-sample adjustment was used} 
#'   \item{cluster}{contains the factor vector that defines independent 
#'   clusters} 
#'   \item{bread}{contains the bread matrix} 
#'   \item{v_scale}{constant used in scaling the sandwich estimator} 
#'   \item{est_mats}{contains
#'   a list of estimating matrices used to calculate the sandwich estimator} 
#'   \item{adjustments}{contains a list of adjustment matrices used to calculate
#'   the sandwich estimator} \item{target}{contains the working
#'   variance-covariance model used to calculate the adjustment matrices. This 
#'   is needed for calculating small-sample corrections for Wald tests.} }
#'   
#' @seealso \code{\link{vcovCR.lm}}, \code{\link{vcovCR.plm}}, 
#'   \code{\link{vcovCR.gls}}, \code{\link{vcovCR.lme}}, 
#'   \code{\link{vcovCR.robu}}, \code{\link{vcovCR.rma.uni}}, 
#'   \code{\link{vcovCR.rma.mv}}
#'   
#' @export
#' @import stats

vcovCR <- function(obj, cluster, type, target, inverse_var, form) UseMethod("vcovCR")

#' Cluster-robust variance-covariance matrix
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @rdname vcovCR
#' @export

vcovCR.default <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich") 
  vcov_CR(obj, cluster, type, target, inverse_var, form)

#---------------------------------------------
# Cluster-robust variance estimator
#---------------------------------------------

adjust_est_mats <- function(type, est_mats, adjustments) {
  switch(type,
         CR0 = est_mats,
         CR1 = lapply(est_mats, function(e) e * adjustments),
         CR1S = lapply(est_mats, function(e) e * adjustments),
         CR2 = Map(function(e, a) e %*% a, e = est_mats, a = adjustments),
         CR3 = Map(function(e, a) e %*% a, e = est_mats, a = adjustments),
         CR4 = Map(function(e, a) a %*% e, e = est_mats, a = adjustments))
}

# uses methods residuals_CS(), model_matrix(), weightMatrix(), 
# targetVariance(), bread(), v_scale()

vcov_CR <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich") {
  
  cluster <- droplevels(as.factor(cluster))
  
  alias <- is.na(coef_CS(obj))
  X <- model_matrix(obj)
  Xp <- projection_matrix(obj)
  if (any(alias)) {
    X <- X[, !alias, drop = FALSE]
    Xp <- Xp[, !alias, drop = FALSE]
  }  
  
  p <- NCOL(X)
  N <- NROW(X)
  
  if (length(cluster) != N) {
    if (class(na.action(obj)) == "omit") {
      cluster <- droplevels(cluster[-na.action(obj)])
    } else {
      stop("Clustering variable must have length equal to nrow(model_matrix(obj)).")
    }
  } 
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  Xp_list <- matrix_list(Xp, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XpW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = Xp_list, w = W_list)
  
  if (is.null(target)) {
    if (inverse_var) {
      Theta_list <- lapply(W_list, function(w) chol2inv(chol(w)))
    } else {
      Theta_list <- targetVariance(obj, cluster)
    }
  } else {
    if (!is.list(target)) {
      Theta_list <- matrix_list(target, cluster, "both")
    } else {
      Theta_list <- target
    }
  }
  
  if (type %in% c("CR2","CR4")) {
    S <- augmented_model_matrix(obj, cluster, inverse_var)
    
    if (is.null(S)) {
      U_list <- Xp_list
      UW_list <- XpW_list
      M_U <- bread(obj) / v_scale(obj)
    } else {
      U <- cbind(Xp, S)
      U_list <- matrix_list(U, cluster, "row")
      UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
      UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
      M_U <- matrix_power(Reduce("+",UWU_list), p = -1)
    }
  }
  
  adjustments <- do.call(type, args = mget(names(formals(type))))
  
  E_list <- adjust_est_mats(type = type, est_mats = XpW_list, adjustments = adjustments)
  
  resid <- residuals_CS(obj)
  res_list <- split(resid, cluster)
  
  components <- do.call(cbind, Map(function(e, r) e %*% r, e = E_list, r = res_list))
  
  v_scale <- v_scale(obj)
  meat <- tcrossprod(components) / v_scale
  
  if (form == "sandwich") {
    bread <- bread(obj)
  } else if (form == "meat") {
    bread <- NULL
  } else if (is.matrix(form)) {
    bread <- form
    form <- "sandwich"
  } 
  
  vcov <- switch(form, 
                 sandwich = bread %*% meat %*% bread / v_scale,
                 meat = meat)
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "type") <- type
  attr(vcov, "cluster") <- cluster
  attr(vcov, "bread") <- bread
  attr(vcov, "v_scale") <- v_scale
  attr(vcov, "est_mats") <- XpW_list
  attr(vcov, "adjustments") <- adjustments
  attr(vcov, "target") <- Theta_list
  attr(vcov, "inverse_var") <- inverse_var
  class(vcov) <- c("vcovCR","clubSandwich")
  return(vcov)
}

#---------------------------------------------
# as.matrix method for vcovCR
#---------------------------------------------

#' @export

as.matrix.clubSandwich <- function(x, ...) {
  attr(x, "type") <- NULL
  attr(x, "cluster") <- NULL
  attr(x, "bread") <- NULL
  attr(x, "v_scale") <- NULL
  attr(x, "est_mats") <- NULL
  attr(x, "adjustments") <- NULL
  attr(x, "target") <- NULL
  attr(x, "inverse_var") <- NULL
  class(x) <- "matrix"
  x
}


#---------------------------------------------
# print method for vcovCR
#---------------------------------------------

#' @export

print.clubSandwich <- function(x, ...) {
  print(as.matrix(x))
}
