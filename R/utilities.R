#----------------------------------------------
# check that CR2 and CR4 are target-unbiased
#----------------------------------------------

check_CR <- function(obj, vcov, ...) {

  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  # calculate E(V^CR)  
  cluster <- attr(vcov, "cluster")
  Theta_list <- attr(vcov, "target")
  S_array <- get_S_array(obj, vcov)
  E_CRj <- lapply(1:nlevels(cluster), function(j) tcrossprod(S_array[,,j]))
         
  # calculate target
  Xp <- projection_matrix(obj)
  alias <- is.na(coef_CR(obj))
  if (any(alias)) Xp <- Xp[, !alias, drop = FALSE]
  p <- NCOL(Xp)
  N <- length(cluster)
  J <- nlevels(cluster)
  
  Xp_list <- matrix_list(Xp, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XpW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = Xp_list, w = W_list)
  XWX_list <- Map(function(xw, x) xw %*% x, xw = XpW_list, x = Xp_list)
  M <- chol2inv(chol(Reduce("+", XWX_list)))
  
  MXWTWXM <- Map(function(xw, theta) M %*% as.matrix(xw %*% theta %*% t(xw)) %*% M, 
                    xw = XpW_list, theta = Theta_list)
  eq <- all.equal(E_CRj, MXWTWXM)
  if (eq==TRUE) TRUE else list(E_CRj = E_CRj, target = MXWTWXM)
}
