#----------------------------------------------
# check that CR2 and CR4 are target-unbiased
#----------------------------------------------

check_CR <- function(obj, vcov, ...) {

  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  # calculate E(V^CR)  
  cluster <- attr(vcov, "cluster")
  E_list <- attr(vcov, "estmats")
  target <- attr(vcov, "target")
  S_array <- get_S_array(obj, cluster, target, E_list)
  E_CRj <- lapply(1:nlevels(cluster), function(j) tcrossprod(S_array[,,j]))
         
  # calculate target
  X <- model_matrix(obj)
  alias <- is.na(coef_CR(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- NCOL(X)
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  W <- weightMatrix(obj)
  W_list <- matrix_list(W, cluster, "both")
  Theta_list <- matrix_list(target, cluster, "both")
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, w = W_list)
  XWX_list <- Map(function(xw, x) xw %*% x, xw = XW_list, x = X_list)
  M <- chol2inv(chol(Reduce("+", XWX_list)))
  
  MXWTWXM <- Map(function(xw, theta) M %*% as.matrix(xw %*% theta %*% t(xw)) %*% M, 
                    xw = XW_list, theta = Theta_list)
  eq <- all.equal(E_CRj, MXWTWXM)
  if (eq==TRUE) TRUE else list(E_CRj = E_CRj, target = MXWTWXM)
}