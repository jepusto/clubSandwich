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
  E_CRj <- apply(S_array, 3, tcrossprod)
  
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
  XW_list <- mapply(function(x, w) as.matrix(t(x) %*% w), 
                     x = X_list, w = W_list, SIMPLIFY = FALSE)
  XW <- matrix(unlist(XW_list), p, N)[,order(order(cluster))]
  M <- chol2inv(chol(XW %*% X))
  MXWTWXM <- mapply(function(x, w, t) M %*% as.matrix(t(x) %*% w %*% t %*% w %*% x) %*% M, 
                    x = X_list, w = W_list, t = Theta_list, SIMPLIFY = TRUE)
  eq <- all.equal(E_CRj, MXWTWXM)
  if (eq==TRUE) TRUE else list(E_CRj = E_CRj, target = MXWTWXM)
}