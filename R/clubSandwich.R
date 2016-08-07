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
#'   
#' @description This is a generic function, with specific methods defined for 
#' \code{\link[stats]{lm}}, \code{\link[plm]{plm}}, \code{\link[nlme]{gls}},
#' \code{\link[nlme]{lme}}, \code{\link[robumeta]{robu}}, \code{\link[metafor]{rma.uni}}, and \code{\link[metafor]{rma.mv}} objects.
#' 
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates. The matrix has several attributes: 
#'   \describe{ \item{type}{indicates which small-sample adjustment was used} 
#'   \item{cluster}{contains the factor vector that defines independent 
#'   clusters} \item{estmats}{contains a list of adjustment matrices used to 
#'   calculate the sandwich estimator, which are needed for calculating 
#'   small-sample corrections for Wald tests} \item{target}{contains the working
#'   variance-covariance model used to calculate the adjustment matrices. This 
#'   is also needed for calculating small-sample corrections for Wald tests.} }
#'   
#' @seealso \code{\link{vcovCR.lm}}, \code{\link{vcovCR.plm}}, 
#'   \code{\link{vcovCR.gls}}, \code{\link{vcovCR.lme}}, 
#'   \code{\link{vcovCR.robu}}, \code{\link{vcovCR.rma.uni}}, 
#'   \code{\link{vcovCR.rma.mv}}
#'   
#' @export
#' @import stats

vcovCR <- function(obj, cluster, type, target, inverse_var) UseMethod("vcovCR")

#' Cluster-robust variance-covariance matrix
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @rdname vcovCR
#' @export

vcovCR.default <- function(obj, cluster, type, target = NULL, inverse_var = FALSE) 
  vcov_CR(obj, cluster, type, target, inverse_var)

#---------------------------------------------
# Cluster-robust variance estimator
#---------------------------------------------

# uses methods residuals_CS(), model_matrix(), weightMatrix(), targetVariance()

vcov_CR <- function(obj, cluster, type, target = NULL, inverse_var = FALSE) {
  
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
  XpWX_list <- Map(function(xw, x) xw %*% x, xw = XpW_list, x = X_list)
  M <- chol2inv(chol(Reduce("+", XpWX_list)))
  
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
      M_U <- M
    } else {
      U <- cbind(Xp, S)
      U_list <- matrix_list(U, cluster, "row")
      UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
      UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
      M_U <- matrix_power(Reduce("+",UWU_list), p = -1)
    }
  }
  
  E_list <- do.call(type, args = mget(names(formals(type))))

  resid <- residuals_CS(obj)

  res_list <- split(resid, cluster)
  
  components <- do.call(cbind, Map(function(e, r) e %*% r, e = E_list, r = res_list))
  vcov <- tcrossprod(components)
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "type") <- type
  attr(vcov, "cluster") <- cluster
  attr(vcov, "estmats") <- E_list
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
  attr(x, "estmats") <- NULL
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

#---------------------------------------------
# matrix manipulation functions
#---------------------------------------------

sub_f <- function(x, fac, dim) {
  function(f) switch(dim,
                      row = x[fac==f, ,drop=FALSE],
                      col = x[ ,fac==f, drop=FALSE],
                      both = x[fac==f, fac==f, drop=FALSE])
}

matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") stop(paste0("Object must be a matrix in order to subset by ",dim,"."))
    x_list <- split(x, fac)
    lapply(x_list, function(x) diag(x, nrow = length(x)))
  } else {
    lapply(levels(fac), sub_f(x, fac, dim)) 
  }
}

matrix_power <- function(x, p, symmetric = TRUE, tol = -12) {
  eig <- eigen(x, symmetric = symmetric)
  val_p <- with(eig, ifelse(values > 10^tol, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

chol_psd <- function(x) with(eigen(x, symmetric=TRUE), sqrt(pmax(values,0)) * t(vectors))

#--------------------------
# get S array
#--------------------------

Sj <- function(e, u, tc, cl, cluster, MUWTheta_cholT) {
  s <- -u %*% MUWTheta_cholT
  s[,cluster==cl] <- tc + s[,cluster==cl]
  e %*% s
}

get_S_array <- function(obj, vcov) {
  
  cluster <- attr(vcov, "cluster")
  E_list <- attr(vcov, "estmats")
  target <- attr(vcov, "target")
  inverse_var <- attr(vcov, "inverse_var")
  
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CS(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- ncol(X)
  
  S <- augmented_model_matrix(obj, cluster, inverse_var)
  
  if (is.null(S)) {
    U <- X
  } else {
    U <- cbind(X, S)
  }
  
  U_list <- matrix_list(U, cluster, "row")
  
  W_list <- weightMatrix(obj, cluster)
  
  UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
  UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
  M_U <- chol2inv(chol(Reduce("+",UWU_list)))

  Theta_cholT <- lapply(target, function(x) t(chol(x)))
  UWThetaC_list <- Map(function(uw, tc) uw %*% tc, uw = UW_list, tc = Theta_cholT)
  MUWTheta_cholT <- M_U %*% (matrix(unlist(UWThetaC_list), ncol(U), N)[,order(order(cluster))])
  
  S_list <- mapply(Sj, e = E_list, u = U_list, tc = Theta_cholT, cl = levels(cluster),
                   MoreArgs = list(cluster=cluster, MUWTheta_cholT=MUWTheta_cholT), SIMPLIFY = FALSE)

  array(unlist(S_list), dim = c(p, N, J))
}
