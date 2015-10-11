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
  
  E_list <- switch(type,
                   CR0 = lapply(XW_list, function(xw) M %*% xw),
                   CR1 = lapply(XW_list, function(xw) (M %*% xw) * sqrt(J / (J - 1))),
                   CR2 = CR2(M, X_list, XW_list, Theta_list = matrix_list(Theta, cluster, "both"), inverse_var),
                   CR3 = CR3(M, X_list, XW_list),
                   CR4 = CR4(M, X_list, XW_list, Theta_list = matrix_list(Theta, cluster, "both"), inverse_var)
                   )

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
# Auxilliary functions for CR* functions
#---------------------------------------------

IH_jj_list <- function(M, X_list, XW_list) {
  mapply(function(x, xw) diag(nrow = nrow(x)) - x %*% M %*% xw,
         x = X_list, xw = XW_list, SIMPLIFY = FALSE)
}

#---------------------------------------------
# Estimating function adjustments
#---------------------------------------------

CR2 <- function(M, X_list, XW_list, Theta_list, inverse_var = FALSE) {

  Theta_chol <- lapply(Theta_list, chol)

  if (inverse_var) {
    IH_jj <- IH_jj_list(M, X_list, XW_list)
    G_list <- mapply(function(a,b,ih) as.matrix(a %*% ih %*% b %*% t(a)), 
                     a = Theta_chol, b = Theta_list, ih = IH_jj, SIMPLIFY = FALSE)
  } else {
    H_jj <- mapply(function(x, xw) x %*% M %*% xw, 
                   x = X_list, xw = XW_list, SIMPLIFY = FALSE)
    xwTwx <- mapply(function(xw, th) xw %*% th %*% t(xw), 
                  xw = XW_list, th = Theta_list, SIMPLIFY = TRUE)
    MXWTWXM <- M %*% matrix(rowSums(xwTwx), nrow(M), ncol(M)) %*% M
    G_list <- mapply(function(thet, h, x, v) 
      as.matrix(v %*% (thet - h %*% thet - thet %*% t(h) + x %*% MXWTWXM %*% t(x)) %*% v),
      thet = Theta_list, h = H_jj, x = X_list, v = Theta_chol, SIMPLIFY = FALSE)
  }

  A_list <- mapply(function(v, g) as.matrix(t(v) %*% Sym_power(g, -1/2) %*% v), 
                   v = Theta_chol, g = G_list, SIMPLIFY = FALSE)

  mapply(function(xw, a) M %*% xw %*% a, xw = XW_list, a = A_list, SIMPLIFY = FALSE)  
}

CR3 <- function(M, X_list, XW_list) {
  IH_jj <- IH_jj_list(M, X_list, XW_list)
  mapply(function(xw, ih) M %*% xw %*% chol2inv(chol(ih)), 
         xw = XW_list, ih = IH_jj, SIMPLIFY = FALSE)
}

CR4 <- function(M, X_list, XW_list, Theta_list, inverse_var = FALSE) {
  
  if (inverse_var) {
    F_list <- mapply(function(xw, x) xw %*% x, 
                     xw = XW_list, x= X_list, SIMPLIFY = FALSE)
    F_chol <- lapply(F_list, chol_psd)
    G_list <- mapply(function(fc, f) fc %*% (f - f %*% M %*% f) %*% t(fc), 
                     fc = F_chol, f = F_list, SIMPLIFY = FALSE)
  } else {
    F_list <- mapply(function(xw, theta) xw %*% theta %*% t(xw), 
                     xw = XW_list, theta = Theta_list, SIMPLIFY = FALSE)
    F_chol <- lapply(F_list, chol_psd)
    XWX_list <- mapply(function(xw, x) xw %*% x, 
                       xw = XW_list, x = X_list, SIMPLIFY = FALSE)
    MXWTWXM <- M %*% apply(array(unlist(F_list), dim = c(dim(M), length(F_list))), 1:2, sum) %*% M
    G_list <- mapply(function(f, fc, xwx)
      as.matrix(fc %*% (f - xwx %*% M %*% f - f %*% M %*% xwx + xwx %*% MXWTWXM %*% xwx) %*% t(fc)),
      f = F_list, fc = F_chol, xwx = XWX_list, SIMPLIFY = FALSE)
  }
  
  D_list <- mapply(function(fc, g) as.matrix(t(fc) %*% Sym_power(g, -1/2) %*% fc), 
                   fc = F_chol, g = G_list, SIMPLIFY = FALSE)
  
  mapply(function(d, xw) M %*% d %*% xw, d = D_list, xw = XW_list, SIMPLIFY = FALSE)
}

#--------------------------
# get S array
#--------------------------

Sj <- function(e, x, tc, cl, cluster, MXWTheta_cholT) {
  s <- -x %*% MXWTheta_cholT
  s[,cluster==cl] <- tc
  e %*% s
}

get_S_array <- function(obj, cluster, target, E_list) {
  
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CR(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- ncol(X)
  X_list <- matrix_list(X, cluster, "row")
  
  W <- weightMatrix(obj)
  W_list <- matrix_list(W, cluster, "both")
  
  XW_list <- mapply(function(x, w) as.matrix(t(x) %*% w), 
                    x = X_list, w = W_list, SIMPLIFY = FALSE)
  XW <- matrix(unlist(XW_list), p, N)[,order(order(cluster))]
  M <- chol2inv(chol(XW %*% X))
  
  if (is.vector(target)) {
    target_sqrt <- sqrt(target)
    Theta_cholT <- matrix_list(target_sqrt, cluster, "both")
    MXWTheta_cholT <- M %*% matrix(unlist(mapply(function(xw, tc) t(tc * t(xw)), 
                                                 xw = XW_list, tc = split(target_sqrt, cluster), SIMPLIFY = FALSE)), 
                                   p, N)
  } else {
    Theta_list <- matrix_list(target, cluster, "both")
    Theta_cholT <- lapply(Theta_list, function(x) t(chol(x)))
    MXWTheta_cholT <- M %*% matrix(unlist(mapply(function(xw, tc) xw %*% tc, xw = XW_list, tc = Theta_cholT, SIMPLIFY = FALSE)), p, N)
  }

  S_list <- mapply(Sj, e = E_list, x = X_list, tc = Theta_cholT, cl = levels(cluster),
                   MoreArgs = list(cluster=cluster, MXWTheta_cholT=MXWTheta_cholT), SIMPLIFY = FALSE)

  array(unlist(S_list), dim = c(p, N, J))
}
