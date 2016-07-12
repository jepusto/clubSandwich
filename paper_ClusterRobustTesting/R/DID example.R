rm(list=ls())
m0 <- 4
m1 <- 14
m_ <- c(m0, m1)
m <- sum(m_)
cluster <- factor(rep(LETTERS[1:m], each = 2))
n <- length(cluster)
time <- rep(c(1,2), m)
trt_clusters <- factor(c(rep(0,m0), rep(1,m1)))
trt <- (time - 1) * rep(trt_clusters, each = 2)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
y <- 0.4 * trt + nu + e

dat <- data.frame(y, time, trt, cluster)
R <- model.matrix(~ 0 + trt, data = dat)
S <- model.matrix(~ 0 + factor(time), data = dat)
T <- model.matrix(~ 0 + cluster, data = dat)

Sp <- residuals(lm.fit(T, S))
Sp2 <- as.matrix(residuals(lm.fit(T, S[,2])))
Rp1 <- as.matrix(residuals(lm.fit(T, R)))
Rp2 <- as.matrix(residuals(lm.fit(Sp, Rp1)))
Rp3 <- as.matrix(residuals(lm.fit(Sp2, Rp1)))
yp1 <- residuals(lm.fit(T, y))
yp2 <- residuals(lm.fit(Sp, yp1))
lm_R <- lm(yp2 ~ 0 + Rp2)
lm_full <- lm(y ~ 0 + cluster + factor(time) + trt, data = dat)

Upf <- residuals(lm.fit(T, cbind(S, R)))
lm_Uf <- lm(yp1 ~ 0 + Upf)

Up1 <- residuals(lm.fit(T, cbind(S[,2], R)))
lm_U1 <- lm(yp1 ~ 0 + Up1)

y_diff <- apply(matrix(y, nrow = 2), 2, diff)
t_Welch <- t.test(y_diff ~ trt_clusters)
d_bar <- tapply(y_diff, trt_clusters, mean)
S_sq <- tapply(y_diff, trt_clusters, var)

all.equal(coef(lm_R)[["Rp2"]], 
          coef(lm_Uf)[["Upftrt"]], 
          coef(lm_U1)[["Up1trt"]], 
          coef(lm_full)[["trt"]])
coef_test(lm_R, vcov = "CR2", cluster = cluster)$SE
coef_test(lm_full, vcov = "CR2", cluster = cluster)["trt","SE"]
coef_test(lm_Uf, vcov = "CR2", cluster = cluster)["Upftrt","SE"]
coef_test(lm_U1, vcov = "CR2", cluster = cluster)["Up1trt","SE"]
sqrt(sum(S_sq / m_))

coef_test(lm_R, vcov = "CR2", cluster = cluster)$df
coef_test(lm_full, vcov = "CR2", cluster = cluster)["trt","df"]
coef_test(lm_Uf, vcov = "CR2", cluster = cluster)["Upftrt","df"]
coef_test(lm_U1, vcov = "CR2", cluster = cluster)["Up1trt","df"]
m^2 * (m0 - 1) * (m1 - 1) / (m0^2 * (m0 - 1) + m1^2 * (m1 - 1))

M_R <- chol2inv(chol(crossprod(Rp2)))
2 * m / (m0 * m1)
M_Uf <- chol2inv(chol(crossprod(Upf)))
M_U1 <- chol2inv(chol(crossprod(Up1)))

Upf_list <- by(Upf, cluster, as.matrix)
Up1_list <- by(Up1, cluster, as.matrix)
R_list <- by(Rp2, cluster, as.matrix)
R_list$A
m1 / (2 * m)
R_list$M
m0 / (2 * m)

IH_R <- IH_jj_list(M = M_R, X_list = R_list, XW_list = lapply(R_list, t))
IH_R$A
1 - m1 / (2 * m0 * m)
IH_R$M
1 - m0 / (2 * m1 * m)

IH_U <- IH_jj_list(M = M_U1, X_list = Up1_list, XW_list = lapply(Up1_list, t))
IH_U$A
1 - 1 / (2 * m0)
IH_U$M
1 - 1 / (2 * m1)

A_mats <- function(X, cluster) {
  M_X <- chol2inv(chol(crossprod(X)))
  X_list <- by(X, cluster, as.matrix)
  IH_jj <- IH_jj_list(M = M_X, X_list = X_list, XW_list = lapply(X_list, t))
  lapply(IH_jj, Sym_power, p = -1/2)
}

A_Upf <- A_mats(X = Upf, cluster)
A_Up1 <- A_mats(X = Up1, cluster)
A_R <- A_mats(X = Rp2, cluster)

mat_neg_sqrt <- function(a) {
  A <- matrix(c(1 - a, a, a, 1 - a), 2, 2)
  A_negsqrt <- Sym_power(A, -1/2)
  list(mat = A_negsqrt, 
       val = (1 + sqrt(1 / (1 - 2 * a))) / 2)
}

mat_neg_sqrt(m1 / (2 * m0 * m))
A_R$A
mat_neg_sqrt(m0 / (2 * m1 * m))
A_R$M

mat_neg_sqrt(1 / (2 * m0))
A_Up1$A
A_Upf$A
mat_neg_sqrt(1 / (2 * m1))
A_Up1$M
A_Upf$M

e_d <- y_diff - d_bar[trt_clusters]
m_i <- m_[trt_clusters]
e_list <- lapply(e_d, function(x) c(-1,1) * x / 2)
V_R <- mapply(function(r, a, e) as.numeric(M_R * t(r) %*% a %*% e)^2, r = R_list, a = A_R, e = e_list)
V_U <- mapply(function(r, a, e) as.numeric(M_R * t(r) %*% a %*% e)^2, r = R_list, a = A_Up1, e = e_list)
V_R * m_i^2 / e_d^2
m_ * m / (m_ * m - rev(m_))

V_U * m_i^2 / e_d^2
m_ / (m_ - 1)

#------------------------------
# Stepped wedge DID
#------------------------------

m <- 12
n <- m + 1
cluster <- factor(rep(LETTERS[1:m], each = n))
N <- length(cluster)
time <- factor(rep(1:n, m))
trt <- as.vector(sapply(1:m, function(x) c(rep(0,x), rep(1, n - x))))
nu <- rnorm(m)[cluster]
e <- rnorm(N)
y <- 0.4 * trt + nu + e

dat <- data.frame(y, time, trt, cluster)
R <- model.matrix(~ 0 + trt, data = dat)
S <- model.matrix(~ 0 + time, data = dat)
T <- model.matrix(~ 0 + cluster, data = dat)
U <- model.matrix(~ 0 + time + trt, data = dat)

Sp <- residuals(lm.fit(T, S))
Up <- residuals(lm.fit(T, U))
Rp1 <- as.matrix(residuals(lm.fit(T, R)))
Rp2 <- as.matrix(residuals(lm.fit(Sp, Rp1)))
yp1 <- residuals(lm.fit(T, y))
yp2 <- residuals(lm.fit(Sp, yp1))
lm_R <- lm(yp2 ~ 0 + Rp2)
lm_U <- lm(yp1 ~ 0 + Up)
lm_full <- lm(y ~ 0 + cluster + time + trt, data = dat)

all.equal(coef(lm_R)[["Rp2"]], 
          coef(lm_U)[["Uptrt"]], 
          coef(lm_full)[["trt"]])
(SE_R <- coef_test(lm_R, vcov = "CR2", cluster = cluster)$SE)
(SE_U <- coef_test(lm_U, vcov = "CR2", cluster = cluster)["Uptrt","SE"])
coef_test(lm_full, vcov = "CR2", cluster = cluster)["trt","SE"]
SE_U / SE_R

(df_R <- coef_test(lm_R, vcov = "CR2", cluster = cluster)$df)
(df_U <- coef_test(lm_U, vcov = "CR2", cluster = cluster)["Uptrt","df"])
coef_test(lm_full, vcov = "CR2", cluster = cluster)["trt","df"]
df_R / df_U
