
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
  ignore_FE <- attr(vcov, "ignore_FE")
  
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CS(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- ncol(X)
  
  W_list <- weightMatrix(obj, cluster)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) w_scale <- 1
  
  S <- augmented_model_matrix(obj, cluster, inverse_var, ignore_FE)
  
  if (is.null(S)) {
    U_list <- matrix_list(X, cluster, "row")
    rm(X, S)
    u <- p
    UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
    M_U <- w_scale * M
  } else {
    U_list <- matrix_list(cbind(X, S), cluster, "row")
    rm(X, S)
    u <- ncol(U_list[[1]])
    UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
    UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
    M_U <- chol2inv(chol(Reduce("+",UWU_list)))
    rm(UWU_list)
  }

  M_U_ct <- t(chol(M_U))
  
  ME_list <- lapply(E_list, function(e) M %*% e)
  
  G_list <- Map(function(me, theta) me %*% t(chol(theta)), me = ME_list, theta = target)
  
  if (inverse_var) {
    H_array <- array(unlist(Map(function(me, u) me %*% u %*% M_U_ct, me = ME_list, u = U_list)), 
                     dim = c(p, u, J))
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
# get P array
#--------------------------

get_P_array <- function(GH, all_terms = FALSE) {
  
  dims <- dim(GH$H)
  
  if (all_terms) {
    
    if (length(dims)==3) {
      P_array <- array(NA, dim = c(dims[1], dims[1], dims[3], dims[3]))
      for (i in 1:dims[1]) for (j in i:dims[1]) {
        if (dims[2] == 1L) {
          tmp <- -tcrossprod(GH$H[i,,], GH$H[j,,])
        } else {
          tmp <- -crossprod(GH$H[i,,], GH$H[j,,])
        }
        diag(tmp) <- diag(tmp) + sapply(GH$G, function(x) sum(x[i,] * x[j,]))
        P_array[i,j,,] <- tmp
        if (j > i) P_array[j,i,,] <- t(tmp)
      }
    } else {
      P_array <- array(NA, dim = c(dims[2], dims[2], dims[4], dims[4]))
      for (i in 1:dims[2]) for (j in i:dims[2]) {
        tmp <- crossprod(GH$H[3,i,,], GH$H[3,j,,]) - 
                  crossprod(GH$H[1,i,,], GH$H[2,j,,]) - 
                    crossprod(GH$H[2,i,,], GH$H[1,j,,])
        diag(tmp) <- diag(tmp) + sapply(GH$G, function(x) sum(x[i,] * x[j,]))
        P_array[i,j,,] <- tmp
        if (j > i) P_array[j,i,,] <- t(tmp)
      }
    }
    
  } else {
    
    if (length(dims)==3) {
      P_array <- array(-apply(GH$H, 1, crossprod), dim = c(dims[3], dims[3], dims[1]))
      P_diag <- matrix(sapply(GH$G, function(x) rowSums(x^2)), nrow = dims[1], ncol = dims[3])
      for (i in 1:dims[1]) diag(P_array[,,i]) <- diag(P_array[,,i]) + P_diag[i,]
    } else {
      P_array <- array(apply(GH$H, 2, function(h) {
        uf <- crossprod(h[1,,], h[2,,])
        crossprod(h[3,,]) - uf - t(uf)
      }), dim = c(dims[4], dims[4], dims[2]))
      P_diag <- matrix(sapply(GH$G, function(x) rowSums(x^2)), nrow = dims[2], ncol = dims[4])
      for (i in 1:dims[2]) diag(P_array[,,i]) <- diag(P_array[,,i]) + P_diag[i,]
    }
    
  }
  
  P_array
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
  ignore_FE <- attr(vcov, "ignore_FE")
  
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X <- model_matrix(obj)
  alias <- is.na(coef_CS(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- ncol(X)
  
  S <- augmented_model_matrix(obj, cluster, inverse_var, ignore_FE)
  
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

