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

vcovCR <- function(obj, cluster, type, target, inverse_var) UseMethod("vcovCR")

vcovCR.default <- function(obj, cluster, type, target = NULL, inverse_var = FALSE) 
  vcov_CR(obj, cluster, type, target, inverse_var)

#---------------------------------------------
# Cluster-robust variance estimator
#---------------------------------------------

# uses methods residuals_CR(), model_matrix(), weightMatrix(), targetVariance()

vcov_CR <- function(obj, cluster, type, target = NULL, inverse_var = FALSE) {
  
  cluster <- droplevels(as.factor(cluster))
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CR(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  
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
  
  resid <- residuals_CR(obj)
  W <- weightMatrix(obj)
  W_list <- matrix_list(W, cluster, "both")
  
  if (is.null(target)) {
    Theta <- targetVariance(obj)
  } else {
    Theta <- target
  }
  
  XW_list <- mapply(function(x, w) as.matrix(t(x) %*% w), 
                    x = X_list, w = W_list, SIMPLIFY = FALSE)
  XW <- matrix(unlist(XW_list), p, N)[,order(order(cluster))]
  M <- chol2inv(chol(XW %*% X))
  
  if (type=="CR0") E_list <- lapply(XW_list, function(x) M %*% x)
  if (type=="CR1") E_list <- lapply(XW_list, function(x) (M %*% x) * J / (J - 1))
  if (type=="CR2") {
    IH <- diag(nrow = N) - X %*% M %*% XW 
    E_list <- CR2(M, XW_list, IH, Theta, cluster, inverse_var)
  }
  if (type=="CR3") {
    IH <- diag(nrow = N) - X %*% M %*% XW
    IH_jj <- mapply(function(x, xw) diag(nrow = nrow(x)) - x %*% M %*% xw,
                    x = X_list, xw = XW_list)
    CR3(M, XW_list, IH_jj = matrix_list(IH, cluster, "both"))
  }
  if (type=="CR4") {
    IH <- diag(nrow = N) - X %*% M %*% XW
    CR4(M, X_list, XW_list, IH, Theta, cluster, inverse_var)
  }
    
  res_list <- split(resid, cluster)
  
  components <- mapply(function(e, r) e %*% r, e = E_list, r = res_list, SIMPLIFY = TRUE)
  
  vcov <- tcrossprod(components)
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "type") <- type
  attr(vcov, "cluster") <- cluster
  attr(vcov, "estmats") <- E_list
  attr(vcov, "target") <- Theta 
  class(vcov) <- c("vcovCR","clubSandwich")
  return(vcov)
}

#---------------------------------------------
# print method for vcovCR
#---------------------------------------------

#' @export

print.clubSandwich <- function(x, ...) {
  attr(x, "type") <- NULL
  attr(x, "cluster") <- NULL
  attr(x, "estmats") <- NULL
  attr(x, "target") <- NULL
  class(x) <- "matrix"
  print(x)
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
    tapply(x, fac, function(x) diag(x, nrow = length(x)))  
  } else {
    lapply(levels(fac), sub_f(x, fac, dim)) 
  }
}

Sym_power <- function(x, p, tol = -12) {
  eig <- eigen(x, symmetric = TRUE)
  val_p <- with(eig, ifelse(values > 10^tol, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

chol_psd <- function(x) with(eigen(x), sqrt(pmax(values,0)) * t(vectors))

#---------------------------------------------
# Estimating function adjustments
#---------------------------------------------

CR2 <- function(M, XW_list, IH, Theta, cluster, inverse_var = FALSE) {

  Theta_list <- matrix_list(Theta, cluster, "both")
  Theta_chol <- lapply(Theta_list, chol)

  if (inverse_var) {
    IH_jj <- matrix_list(IH, cluster, "both")
    G_list <- mapply(function(a,b,ih) as.matrix(a %*% ih %*% b %*% t(a)), 
                     a = Theta_chol, b = Theta_list, ih = IH_jj, SIMPLIFY = FALSE)
  } else {
    IH_list <- matrix_list(IH, cluster, "row")
    if (is.vector(Theta)) {
      G_list <- mapply(function(v, ih) as.matrix(v %*% ih %*% (Theta * t(ih)) %*% t(v)),
                       v = Theta_chol, ih = IH_list, SIMPLIFY = FALSE)
    } else {
      G_list <- mapply(function(v, ih) as.matrix(v %*% ih %*% Theta %*% t(ih) %*% t(v)), 
                       v = Theta_chol, ih = IH_list, SIMPLIFY = FALSE)
    }
  }

  A_list <- mapply(function(v, g) as.matrix(t(v) %*% Sym_power(g, -1/2) %*% v), 
                   v = Theta_chol, g = G_list, SIMPLIFY = FALSE)

  mapply(function(xw, a) M %*% xw %*% a, xw = XW_list, a = A_list, SIMPLIFY = FALSE)  
}

CR3 <- function(M, XW_list, IH_jj) {
  mapply(function(xw, ih) M %*% xw %*% chol2inv(chol(ih)), 
         xw = XW_list, ih = IH_jj, SIMPLIFY = FALSE)
}

CR4 <- function(M, X_list, XW_list, IH, Theta, cluster, inverse_var = FALSE) {
  
  if (inverse_var) {
    F_list <- mapply(function(xw, x) xw %*% x, 
                     xw = XW_list, x= X_list, SIMPLIFY = FALSE)
    F_chol <- lapply(F_list, chol_psd)
    G_list <- mapply(function(f_c, f) f_c %*% (f - f %*% M %*% f) %*% t(f_c), 
                     f_c = F_chol, f = F_list, SIMPLIFY = FALSE)
  } else {
    Theta_list <- matrix_list(Theta, cluster, "both")
    F_list <- mapply(function(xw, theta) xw %*% theta %*% t(xw), 
                     xw = XW_list, theta = Theta_list, SIMPLIFY = FALSE)
    F_chol <- lapply(F_list, chol_psd)
    IH_list <- matrix_list(IH, cluster, "row")
    XWIH_list <- mapply(function(xw, ih) as.matrix(xw %*% ih), 
                        xw = XW_list, ih = IH_list, SIMPLIFY = FALSE)
    if (is.vector(Theta)) {
      G_list <- mapply(function(f, xwih) as.matrix(f %*% xwih %*% (Theta * t(xwih)) %*% t(f)), 
                       f = F_chol, xwih = XWIH_list, SIMPLIFY = FALSE)
    } else {
      G_list <- mapply(function(f, xwih) as.matrix(f %*% xwih %*% Theta %*% t(xwih) %*% t(f)), 
                       f = F_chol, xwih = XWIH_list, SIMPLIFY = FALSE)
    }
  }
  
  D_list <- mapply(function(f, g) as.matrix(t(f) %*% Sym_power(g, -1/2) %*% f), 
                   f = F_chol, g = G_list, SIMPLIFY = FALSE)
  
  mapply(function(d, xw) M %*% d %*% xw, d = D_list, xw = XW_list, SIMPLIFY = FALSE)
}


#--------------------------
# get S array
#--------------------------

get_S_array <- function(obj, cluster, target, E_list) {
  
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CR(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  W <- weightMatrix(obj)
  if (is.vector(W)) {
    XW <- t(W * X)
  } else {
    XW <- t(X) %*% W
  }
  M <- chol2inv(chol(XW %*% X))
  
  IH <- as.matrix(diag(nrow = NROW(X)) - X %*% M %*% XW)
  IH_list <- matrix_list(IH, cluster, "row")
  
  if (is.vector(target)) {
    Theta_cholT <- sqrt(target)
    S_list <- mapply(function(e, ih) e %*% (ih * matrix(Theta_cholT, nrow = nrow(ih), ncol = ncol(ih), byrow = TRUE)),
                     e = E_list, ih = IH_list, SIMPLIFY = FALSE)
  } else {
    Theta_cholT <- t(chol(target))
    S_list <- mapply(function(e, ih) as.matrix(e %*% ih %*% Theta_cholT),
                     e = E_list, ih = IH_list, SIMPLIFY = FALSE)
  }
  
  array(unlist(S_list), dim = c(ncol(X), N, J))
}
