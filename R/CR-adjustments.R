#---------------------------------------------
# Auxilliary functions for CR* functions
#---------------------------------------------

IH_jj_list <- function(M, X_list, XW_list) {
  Map(function(x, xw) diag(nrow = nrow(x)) - x %*% M %*% xw,
         x = X_list, xw = XW_list)
}

#---------------------------------------------
# Estimating function adjustments
#---------------------------------------------

CR0 <- function(XpW_list, M) 
  lapply(XpW_list, function(xw) M %*% xw)

CR1 <- function(XpW_list, M, J) 
  lapply(XpW_list, function(xw) (M %*% xw) * sqrt(J / (J - 1)))

CR1S <- function(XpW_list, M, J, N, p) 
  lapply(XpW_list, function(xw) (M %*% xw) * sqrt(J * N / ((J - 1) * (N - p))))

CR2 <- function(M_U, U_list, UW_list, M, XpW_list, X_list, Theta_list, inverse_var = FALSE) {
  
  Theta_chol <- lapply(Theta_list, chol)
  
  if (inverse_var) {
    IH_jj <- IH_jj_list(M_U, U_list, UW_list)
    G_list <- Map(function(a,b,ih) as.matrix(a %*% ih %*% b %*% t(a)), 
                     a = Theta_chol, b = Theta_list, ih = IH_jj)
  } else {
    H_jj <- Map(function(u, uw) u %*% M_U %*% uw, 
                   u = U_list, uw = UW_list)
    uwTwu <- Map(function(uw, th) uw %*% th %*% t(uw), 
                    uw = UW_list, th = Theta_list)
    MUWTWUM <- M_U %*% Reduce("+", uwTwu) %*% M_U
    G_list <- Map(function(thet, h, u, v) 
      as.matrix(v %*% (thet - h %*% thet - thet %*% t(h) + u %*% MUWTWUM %*% t(u)) %*% t(v)),
      thet = Theta_list, h = H_jj, u = U_list, v = Theta_chol)
  }
  
  A_list <- Map(function(v, g) as.matrix(t(v) %*% Sym_power(g, -1/2) %*% v), 
                   v = Theta_chol, g = G_list)
  
  Map(function(xw, a) M %*% xw %*% a, xw = XpW_list, a = A_list)  
}

CR3 <- function(M, X_list, XpW_list) {
  IH_jj <- IH_jj_list(M, X_list, XpW_list)
  Map(function(xw, ih) M %*% xw %*% solve(ih), 
         xw = XpW_list, ih = IH_jj)
}

CR4 <- function(M, X_list, XpW_list, Theta_list, inverse_var = FALSE) {
  
  if (inverse_var) {
    F_list <- Map(function(xw, x) xw %*% x, xw = XpW_list, x= X_list)
    F_chol <- lapply(F_list, chol_psd)
    G_list <- Map(function(fc, fm) fc %*% (fm - fm %*% M %*% fm) %*% t(fc), 
                  fc = F_chol, fm = F_list)
  } else {
    F_list <- Map(function(xw, theta) xw %*% theta %*% t(xw), 
                     xw = XpW_list, theta = Theta_list)
    F_chol <- lapply(F_list, chol_psd)
    XWX_list <- Map(function(xw, x) xw %*% x, xw = XpW_list, x = X_list)
    MXWTWXM <- M %*% Reduce("+", F_list) %*% M
    G_list <- Map(function(fm, fc, xwx)
      as.matrix(fc %*% (fm - xwx %*% M %*% fm - fm %*% M %*% xwx + xwx %*% MXWTWXM %*% xwx) %*% t(fc)),
      fm = F_list, fc = F_chol, xwx = XWX_list)
  }
  
  D_list <- Map(function(fc, g) as.matrix(t(fc) %*% Sym_power(g, -1/2) %*% fc), 
                   fc = F_chol, g = G_list)
  
  Map(function(d, xw) M %*% d %*% xw, d = D_list, xw = XpW_list)
}
