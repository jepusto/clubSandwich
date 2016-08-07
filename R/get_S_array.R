#--------------------------
# get S array
#--------------------------

Sj <- function(M, e, u, tc, cl, cluster, MUWTheta_cholT) {
  s <- -u %*% MUWTheta_cholT
  s[,cluster==cl] <- tc + s[,cluster==cl]
  M %*% e %*% s
}

get_S_array <- function(obj, vcov) {
  
  cluster <- attr(vcov, "cluster")
  M <- attr(vcov, "bread") / attr(vcov, "v_scale")

  E_list <- adjust_est_mats(type = attr(vcov, "type"), 
                            est_mats = attr(vcov, "est_mats"), 
                            adjustments = attr(vcov, "adjustments"))
  
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
                   MoreArgs = list(M = M, cluster=cluster, MUWTheta_cholT=MUWTheta_cholT), 
                   SIMPLIFY = FALSE)
  
  array(unlist(S_list), dim = c(p, N, J))
}
