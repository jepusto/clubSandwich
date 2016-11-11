#--------------------------
# get P array
#--------------------------

get_P_array <- function(obj, vcov, verbose = FALSE) {
  
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
  rm(S)
  U_list <- matrix_list(U, cluster, "row")
  rm(U)
  
  W_list <- weightMatrix(obj, cluster)
  
  UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
  UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
  M_U <- chol2inv(chol(Reduce("+",UWU_list)))
  rm(UWU_list)
  
  EU_list <- Map(function(e, u) e %*% u, e = E_list, u = U_list)
  
  P_array <- array(NA, dim = c(p, p, J, J))
  
  if (verbose) pb <- txtProgressBar(min = 0, max = J^2, style = 3)
  
  if (inverse_var) {
    for (i in 1:J) {
      P_array[,,i,i] <- E_list[[i]] %*% target[[i]] %*% t(E_list[[i]]) - EU_list[[i]] %*% M_U %*% t(EU_list[[i]])
      if (i < J) {
        for (j in (i + 1):J) {
          if (verbose) setTxtProgressBar(pb, J * (i - 1) + j)
          p <- -EU_list[[i]] %*% M_U %*% t(EU_list[[j]])
          P_array[,,i,j] <- p
          P_array[,,j,i] <- t(p)
        }
      }
    }
  } else {
    TWU_list <- Map(function(t, w, u) t %*% w %*% u, t = target, w = W_list, u = U_list)
    EF_list <- Map(function(e, twu) e %*% twu, e = E_list, twu = TWU_list)
    UWTWU_list <- Map(function(uw, twu) uw %*% twu, uw = UW_list, twu = TWU_list)
    Omega <- M_U %*% Reduce("+", UWTWU_list) %*% M_U
    rm(TWU_list, UWTWU_list)
    
    for (i in 1:J) {
      P_array[,,i,i] <- E_list[[i]] %*% target[[i]] %*% t(E_list[[i]]) + EU_list[[i]] %*% Omega %*% t(EU_list[[i]]) - EU_list[[i]] %*% M_U %*% t(EF_list[[i]]) - EF_list[[i]] %*% M_U %*% t(EU_list[[i]])
      if (i < J) {
        for (j in (i + 1):J) {
          if (verbose) setTxtProgressBar(pb, J * (i - 1) + j)
          p <- EU_list[[i]] %*% Omega %*% t(EU_list[[j]]) - EU_list[[i]] %*% M_U %*% t(EF_list[[j]]) - EF_list[[i]] %*% M_U %*% t(EU_list[[j]])
          P_array[,,i,j] <- p
          P_array[,,j,i] <- t(p)
        }
      }
    }
  }
  
  array(apply(P_array, 3:4, function(x) M %*% x %*% M), dim = dim(P_array))
}

#--------------------------
# get G list and H array
#--------------------------

get_GH <- function(obj, vcov) {
  
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
  rm(S)
  U_list <- matrix_list(U, cluster, "row")
  u <- ncol(U)
  rm(U)
  
  W_list <- weightMatrix(obj, cluster)
  
  UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
  UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
  M_U <- chol2inv(chol(Reduce("+",UWU_list)))
  M_U_ct <- t(chol(M_U))
  rm(UWU_list)
  
  ME_list <- lapply(E_list, function(e) M %*% e)
  
  G_list <- Map(function(me, theta) me %*% t(chol(theta)), me = ME_list, theta = target)
  
  
  if (inverse_var) {
    H_array <- array(Map(function(me, u) me %*% u %*% M_U_ct, me = ME_list, u = U_list), dim = c(p, u, J))
  } else {
    H_array <- array(NA, dim = c(3, p, u, J))
    MEU_list <- Map(function(me, u) me %*% u, me = ME_list, u = U_list)
    H_array[1,,,] <- unlist(lapply(MEU_list, function(meu) meu %*% M_U_ct))
    TWU_list <- Map(function(t, w, u) t %*% w %*% u, t = target, w = W_list, u = U_list)
    MEF_list <- Map(function(me, twu) me %*% twu, me = ME_list, twu = TWU_list)
    H_array[2,,,] <- unlist(lapply(MEF_list, function(mef) mef %*% M_U_ct))
    rm(MEF_list)
    UWTWU_list <- Map(function(uw, twu) uw %*% twu, uw = UW_list, twu = TWU_list)
    Omega_ct <- t(chol(M_U %*% Reduce("+", UWTWU_list) %*% M_U))
    rm(TWU_list, UWTWU_list)
    H_array[3,,,] <- unlist(lapply(MEU_list, function(meu) meu %*% Omega_ct))
    rm(MEU_list, Omega_ct)
  }
  
  list(G = G_list, H = H_array)
}

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

