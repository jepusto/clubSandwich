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

CR0 <- function(XW_list, M) 
  lapply(XW_list, function(xw) M %*% xw)

CR1 <- function(XW_list, M, J) 
  lapply(XW_list, function(xw) (M %*% xw) * sqrt(J / (J - 1)))

CR1S <- function(XW_list, M, J, N, p) 
  lapply(XW_list, function(xw) (M %*% xw) * sqrt(J * N / ((J - 1) * (N - p))))

CR2 <- function(M_U, U_list, UW_list, M, XW_list, Theta_list, inverse_var = FALSE) {
  
  Theta_chol <- lapply(Theta_list, chol)
  
  if (inverse_var) {
    IH_jj <- IH_jj_list(M_U, U_list, UW_list)
    G_list <- mapply(function(a,b,ih) as.matrix(a %*% ih %*% b %*% t(a)), 
                     a = Theta_chol, b = Theta_list, ih = IH_jj, SIMPLIFY = FALSE)
  } else {
    H_jj <- mapply(function(u, uw) u %*% M_U %*% uw, 
                   u = U_list, uw = UW_list, SIMPLIFY = FALSE)
    uwTwu <- mapply(function(uw, th) uw %*% th %*% t(uw), 
                    uw = UW_list, th = Theta_list, SIMPLIFY = TRUE)
    MUWTWUM <- M_U %*% matrix(rowSums(uwTwu), nrow(M), ncol(M)) %*% M_U
    G_list <- mapply(function(thet, h, u, v) 
      as.matrix(v %*% (thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u)) %*% v),
      thet = Theta_list, h = H_jj, u = U_list, v = Theta_chol, SIMPLIFY = FALSE)
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
