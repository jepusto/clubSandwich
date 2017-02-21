#-----------------------------------------------------
# check that bread can be re-constructed from X and W
#-----------------------------------------------------

check_bread <- function(obj, cluster, y, tol = .Machine$double.eps^0.5) {
  cluster <- droplevels(as.factor(cluster))
  B <- sandwich::bread(obj) / v_scale(obj)
  X_list <- matrix_list(model_matrix(obj), cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XWX <- Reduce("+", Map(function(x, w) t(x) %*% w %*% x, x = X_list, w = W_list))
  M <- chol2inv(chol(XWX))
  attr(M, "dimnames") <- attr(B, "dimnames")
  
  coef <- coef_CS(obj)
  y_list <- split(y, cluster)
  XWy <- Reduce("+", Map(function(x, w, y) t(x) %*% w %*% y, x = X_list, w = W_list, y = y_list))
  beta <- as.vector(M %*% XWy)
  names(beta) <- names(coef)
  
  eq_bread <- diff(range(M / B)) < tol
  eq_coef <- all.equal(beta, coef, tol = tol)
  if (all(c(eq_coef, eq_bread) == TRUE)) TRUE else list(M = M, B = B, beta = beta, coef = coef)
}

#----------------------------------------------
# check that CR2 and CR4 are target-unbiased
#----------------------------------------------

check_CR <- function(obj, vcov, ..., tol = .Machine$double.eps^0.5) {

  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  # calculate E(V^CRj)  
  cluster <- attr(vcov, "cluster")
  S_array <- get_S_array(obj, vcov)
  E_CRj <- lapply(1:nlevels(cluster), function(j) tcrossprod(S_array[,,j]))
         
  # calculate target
  Theta_list <- attr(vcov, "target")
  Xp <- projection_matrix(obj)
  alias <- is.na(coef_CS(obj))
  if (any(alias)) Xp <- Xp[, !alias, drop = FALSE]
  p <- NCOL(Xp)
  N <- length(cluster)
  J <- nlevels(cluster)
  
  Xp_list <- matrix_list(Xp, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XpW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = Xp_list, w = W_list)
  M <- attr(vcov, "bread") / attr(vcov, "v_scale")
  attr(M, "dimnames") <- NULL

  MXWTWXM <- Map(function(xw, theta) M %*% as.matrix(xw %*% theta %*% t(xw)) %*% M, 
                    xw = XpW_list, theta = Theta_list)
  eq <- all.equal(E_CRj, MXWTWXM, tolerance = tol)
  if (all(eq==TRUE)) TRUE else list(E_CRj = E_CRj, target = MXWTWXM)
}
