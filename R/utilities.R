#----------------------------------------------
# check that CR2 and CR4 are target-unbiased
#----------------------------------------------

check_CR <- function(obj, cluster, target = NULL, inverse_var = FALSE) {
  
  cluster <- as.factor(cluster)
  
  X <- model.matrix(obj)
  alias <- is.na(coef(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  X_list <- matrix_list(X, cluster, "row")
  
  resid <- residuals(obj)
  W <- weightMatrix(obj)
  W_list <- matrix_list(W, cluster, "both")
  
  if (is.null(target)) {
    Theta <- targetVariance(obj)
  } else {
    Theta <- target
  }
  
  J <- nlevels(cluster)
  p <- NCOL(X)
  N <- NROW(X)
  
  XW_list <- mapply(function(x, w) as.matrix(t(x) %*% w), 
                    x = X_list, w = W_list, SIMPLIFY = FALSE)
  XW <- matrix(unlist(XW_list), p, N)[,order(order(cluster))]
  M <- chol2inv(chol(XW %*% X))
  IH <- diag(nrow = N) - X %*% M %*% XW 
  
  # get adjustment functions
  E_CR2 <- CR2(M = diag(1L, nrow = p), XW_list, IH, Theta, cluster, inverse_var)
  E_CR4 <- CR4(M = diag(1L, nrow = p), X_list, XW_list, IH, Theta, cluster, inverse_var)
  
  # get targets
  Theta_list <- matrix_list(Theta, cluster, "both")
  F_list <- mapply(function(xw, theta) as.matrix(xw %*% theta %*% t(xw)), xw = XW_list, theta = Theta_list, SIMPLIFY = FALSE)
  
  # check cluster-by-cluster adjustments  
  Theta_cholT <- t(chol(Theta))
  IH_list <- matrix_list(IH, cluster, "row")
  S_CR2 <- mapply(function(e, ih) as.matrix(e %*% ih %*% Theta_cholT), e = E_CR2, ih = IH_list, SIMPLIFY = FALSE)
  shot_CR2 <- lapply(S_CR2, tcrossprod)
  check_CR2 <- mapply(function(target, current) all.equal(target, current, check.attributes = FALSE), 
                      target = F_list, current = shot_CR2)
  
  S_CR4 <- mapply(function(e, ih) as.matrix(e %*% ih %*% Theta_cholT), e = E_CR4, ih = IH_list, SIMPLIFY = FALSE)
  shot_CR4 <- lapply(S_CR4, tcrossprod)
  check_CR4 <- mapply(function(target, current) all.equal(target, current, check.attributes = FALSE), 
                      target = F_list, current = shot_CR4)
  
  checks <- data.frame(cluster = levels(cluster), CR2 = check_CR2, CR4 = check_CR4)
  
  results <- list(checks = checks, targets = F_list, shot_CR2 = shot_CR2, shot_CR4 = shot_CR4)
  class(results) <- "check_clubSandwich"
  return(results)
}

print.check_clubSandwich <- function(x, ...) {
  print(x$checks, row.names = FALSE)
}

#-------------------------------------------
# Get CR2 and CR4 adjustment matrices
#-------------------------------------------

CR2_A <- function(obj, cluster, target = NULL, inverse_var = FALSE) {
  
  cluster <- droplevels(as.factor(cluster))
  
  X <- model.matrix(obj)
  alias <- is.na(coef(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  
  p <- NCOL(X)
  N <- NROW(X)
  
  if (length(cluster) != N) {
    if (class(na.action(obj)) == "omit") {
      cluster <- droplevels(cluster[-na.action(obj)])
    } else {
      stop("Clustering variable must have length equal to nrow(model.matrix(obj)).")
    }
  } 
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  
  resid <- residuals(obj)
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
  IH <- diag(nrow = N) - X %*% M %*% XW 
  
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
  
  mapply(function(v, g) as.matrix(t(v) %*% Sym_power(g, -1/2) %*% v), 
                   v = Theta_chol, g = G_list, SIMPLIFY = FALSE)
}


CR4_A <- function(obj, cluster, target = NULL, inverse_var = FALSE) {
  
  cluster <- droplevels(as.factor(cluster))
  
  X <- model.matrix(obj)
  alias <- is.na(coef(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  
  p <- NCOL(X)
  N <- NROW(X)
  
  if (length(cluster) != N) {
    if (class(na.action(obj)) == "omit") {
      cluster <- droplevels(cluster[-na.action(obj)])
    } else {
      stop("Clustering variable must have length equal to nrow(model.matrix(obj)).")
    }
  } 
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  
  resid <- residuals(obj)
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
  IH <- diag(nrow = N) - X %*% M %*% XW 
  
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
  
  mapply(function(f, g) as.matrix(t(f) %*% Sym_power(g, -1/2) %*% f), 
                   f = F_chol, g = G_list, SIMPLIFY = FALSE)
}