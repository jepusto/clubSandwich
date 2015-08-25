library(plm)
rm(list=ls())
devtools::load_all()
set.seed(12)

m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
table(cluster)
n <- length(cluster)
X <- matrix(rnorm(3 * n), n, 3)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
w <- rgamma(n, shape = 3, scale = 3)
y <- X %*% c(.4, .3, -.3) + nu + e

dat <- data.frame(y, X, cluster)

lm_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat)
plm_fit <- plm(y ~ X1 + X2 + X3, data = dat, effect = "individual", model = "within", index = c("cluster"))

WLS_fit <- lm(y ~ 0 + cluster + X1 + X2 + X3, data = dat, weights = w)
S <- model.matrix(~ 0 + cluster, data = dat)
M_S <- chol2inv(chol(t(S) %*% (w * S)))
IH_S <- diag(n) - S %*% M_S %*% t(w * S)
y_ <- IH_S %*% y
R_ <- IH_S %*% X
dat_ <- data.frame(y = y_, R_, cluster)
WLS_absorb <- lm(y ~ 0 + X1 + X2 + X3, data = dat_, weights = w)

# coefficient estimates are identical
cbind(coef(lm_fit)[m + 1:3], coef(plm_fit))
all.equal(coef(lm_fit)[m + 1:3], coef(plm_fit))

cbind(coef(WLS_fit)[m + 1:3], coef(WLS_absorb))
all.equal(coef(WLS_fit)[m + 1:3], coef(WLS_absorb))

# CR0 - unweighted
CR0_lm <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR0")
CR0_plm <- vcovCR(plm_fit, type = "CR0")
CR0_lm[m + 1:3, m + 1:3]
CR0_plm
all.equal(CR0_lm[m + 1:3, m + 1:3], CR0_plm[1:3,1:3])

# CR0 - weighted
CR0_WLS <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR0")
CR0_absorb <- vcovCR(WLS_absorb, cluster = dat$cluster, type = "CR0")
CR0_WLS[m + 1:3, m + 1:3]
CR0_absorb
all.equal(CR0_WLS[m + 1:3, m + 1:3], CR0_absorb[1:3,1:3])

# CR2 - unweighted
CR2_lm <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
all.equal(CR2_lm[m + 1:3,m + 1:3], vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = FALSE)[m + 1:3,m + 1:3])
CR2_plm <- vcovCR(plm_fit, type = "CR2")
all.equal(CR2_plm, vcovCR(plm_fit, type = "CR2", inverse_var = FALSE))
CR2_lm[m + 1:3, m + 1:3]
CR2_plm
all.equal(CR2_lm[m + 1:3, m + 1:3], CR2_plm[1:3,1:3])

coef_test(lm_fit, vcov = CR2_lm)[m + 1:3,]
coef_test(plm_fit, vcov = CR2_plm)

Wald_test(lm_fit, constraints = m + 1:3, vcov = CR2_lm, test = "All")
Wald_test(plm_fit, constraints = 1:3, vcov = CR2_plm, test = "All")

# CR2 - weighted
CR2_WLS <- vcovCR(WLS_fit, cluster = dat$cluster, type = "CR2", target = 1 / w, inverse_var = TRUE)
CR2_absorb <- vcovCR(WLS_absorb, cluster = dat$cluster, type = "CR2", target = 1 / w, inverse_var = TRUE)

CR2_WLS[m + 1:3, m + 1:3]
CR2_absorb
all.equal(CR2_WLS[m + 1:3, m + 1:3], CR2_absorb[1:3,1:3]) # differences due to numerical imprecision

coef_test(WLS_fit, vcov = CR2_WLS)[m + 1:3,]
coef_test(WLS_absorb, vcov = CR2_absorb)

Wald_test(WLS_fit, constraints = m + 1:3, vcov = CR2_WLS, test = "All")
Wald_test(WLS_absorb, constraints = 1:3, vcov = CR2_absorb, test = "All")


#------------------------------
# checking my algebra
#------------------------------

X <- model.matrix(~ 0 + cluster + X1 + X2 + X3, data = dat)
M_X <- chol2inv(chol(t(X) %*% (w * X)))
IH_X <- diag(n) - X %*% M_X %*% t(w * X)
M_R <- chol2inv(chol(t(R_) %*% (w * R_)))
IH_R <- diag(n) - R_ %*% M_R %*% t(w * R_) 
all.equal(IH_X, IH_R %*% IH_S)

D <- IH_R %*% IH_S %*% diag(1 / w) %*% t(IH_S) %*% t(IH_R)
Dj_list <- lapply(levels(cluster), function(c) D[cluster==c, cluster==c])
Dj_check <- lapply(levels(cluster), function(c)
                   diag(1 / w[cluster==c]) 
                   - R_[cluster==c,] %*% M_R %*% t(R_[cluster==c,]) 
                   - S[cluster==c,] %*% M_S %*% t(S[cluster==c,]))
mapply(all.equal, Dj_list, Dj_check)
Dj_inv <- lapply(Dj_list, Sym_power, p = -1)

Sj <- lapply(levels(cluster), function(c) S[cluster==c,])
Rj <- lapply(levels(cluster), function(c) R_[cluster==c,])
wj <- lapply(levels(cluster), function(c) w[cluster==c])
RWS <- mapply(function(r, w, s) t(r) %*% (w * s), r = Rj, w = wj, s = Sj, SIMPLIFY = FALSE)
lapply(RWS, round, 12)

U_j <- mapply(function(r, wi) chol2inv(chol(diag(1 / wi) - r %*% M_R %*% t(r))), r = Rj, wi = wj, SIMPLIFY = FALSE)

Uj_Sj <- mapply(function(u, s) u %*% s, u = U_j, s = Sj)
wj_Sj <- mapply(function(wi, s) wi * s, wi = wj, s = Sj)
mapply(all.equal, Uj_Sj, wj_Sj, check.attributes = FALSE)

j <- 1
eig_j <- eigen(M_S - M_S %*% t(Sj[[j]]) %*% (wj[[j]] * Sj[[j]]) %*% M_S, symmetric = TRUE)
with(eig_j, vectors %*% (ifelse(values > 10^-12, values^-1, 0) * t(vectors)))
Sym_power(M_S - M_S %*% t(Sj[[j]]) %*% (wj[[j]] * Sj[[j]]) %*% M_S, p = -1)

(SBZ_eig <- mapply(function(s, wi) eigen(M_S - M_S %*% t(s) %*% (wi * s) %*% M_S, symmetric=TRUE), s = Sj, wi = wj, SIMPLIFY = FALSE) )
SBZ_j <- mapply(function(s, wi) (wi * s) %*% M_S %*% Sym_power(M_S - M_S %*% t(s) %*% (wi * s) %*% M_S, -1) %*% M_S %*% t(wi * s), 
                s = Sj, wi = wj, SIMPLIFY = FALSE)
SBZ_j
E <- IH_R %*% diag(1 / w) %*% t(IH_R)
Ej_inv <- lapply(levels(cluster), function(c) Sym_power(E[cluster==c, cluster==c], -1))
mapply(all.equal, U_j, Ej_inv)
mapply(all.equal, U_j, Dj_inv)
