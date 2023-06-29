set.seed(20230114)

ni <- c(2,3,5)
m <- length(ni)
ri <- lapply(ni, \(x) 1:x)
wi <- lapply(ni, \(x) 1 / (1:x))
rddi_ols <- lapply(ri, \(r) r - mean(r))
rddi_wt <- mapply(\(r, w) r - weighted.mean(r,w = w), r = ri, w = wi)

clust <- factor(rep(LETTERS[1:m], ni))
mu <- 0.2 * unlist(ri) + ni[clust]
yi <- round(rnorm(sum(ni), mean = mu, sd = sqrt(unlist(ri))), 1)
yddi_ols <- tapply(yi, clust, \(x) x - mean(x))
yddi_wt <- mapply(\(y,w) y - weighted.mean(y, w = w), y = split(yi, clust), w = wi)

dat <- data.frame(
  y = yi, 
  y_dd_ols = unlist(yddi_ols),
  y_dd_wt = unlist(yddi_wt),
  R = unlist(ri),
  R_dd_ols = unlist(rddi_ols),
  R_dd_wt = unlist(rddi_wt),
  w = unlist(wi),
  clust = clust
)

matrix_power <- function(x, p, symmetric = TRUE, tol = -12) {
  eig <- eigen(x, symmetric = symmetric)
  val_p <- with(eig, ifelse(values > 10^tol, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

compute_AWR <- function(mod, cluster, target = NULL, R = model.matrix(mod)) {
  
  X <- model.matrix(mod)
  w <- weights(mod)
  if (is.null(w)) w <- rep(1, nrow(X))
  if (is.null(target)) target <- 1 / w

  wX <- w * X
  M <- solve(crossprod(X, wX))
  MXW <- M %*% t(wX)
  I_H <- diag(rep(1, nrow(X))) - X %*% MXW
  I_H_i <- by(I_H, cluster, as.matrix)
  D_i <- split(sqrt(target), cluster)
  
  B_i <- mapply(\(ih, d) tcrossprod(d) * (ih %*% (target * t(ih))),
                ih = I_H_i, d = D_i)
  A_i <- mapply(\(b, d) tcrossprod(d) * matrix_power(b, p = -1/2), 
                b = B_i, d = D_i)
  
  WR_i <- by(w * R, cluster, as.matrix)
  
  mapply(\(a, wr) a %*% wr, a = A_i, wr = WR_i)
}

compute_VCR <- function(mod, cluster, target = NULL, R = model.matrix(mod), include_AWR = TRUE) {
  resids <- residuals(mod)
  w <- weights(mod)
  if (is.null(w)) w <- rep(1, length(resids))
  e_i <- split(resids, cluster)
  
  MR <- solve(crossprod(R, w * R))
  
  AWR_i <- compute_AWR(mod, cluster, target = target, R = R)
  eAWR_i <- mapply(crossprod, x = e_i, y = AWR_i, SIMPLIFY = FALSE)
  meat <- tcrossprod(do.call(cbind, eAWR_i))
  
  V_CR <- MR %*% meat %*% MR
  
  if (include_AWR) {
    rbind(do.call(rbind, AWR_i), V_CR)
  } else{
    V_CR
  }
}

# Ordinary least squares, ID working model

ols_fit <- lm(y ~ clust + R, data = dat)
res_ols_full <- compute_VCR(mod = ols_fit, cluster = dat$clust, R = dat$R_dd_ols)

ols_absorb <- lm(y_dd_ols ~ 0 + R_dd_ols, data = dat)
res_ols_absorb <- compute_VCR(mod = ols_absorb, cluster = dat$clust, R = dat$R_dd_ols)


# Ordinary least squares, heteroskedastic working model

res_ols_full_het <- compute_VCR(mod = ols_fit, cluster = dat$clust, target = dat$R, R = dat$R_dd_ols)
res_ols_absorb_het <- compute_VCR(mod = ols_absorb, cluster = dat$clust, target = dat$R, R = dat$R_dd_ols)

# Weighted least squares, heteroskedastic working model

wls_fit <- lm(y ~ clust + R, weights = w, data = dat)
res_wls_full <- compute_VCR(mod = wls_fit, cluster = dat$clust, R = dat$R_dd_wt)

wls_absorbed <- lm(y_dd_wt ~ 0 + R_dd_wt, weights = w, data = dat)
res_wls_absorb <- compute_VCR(mod = wls_absorbed, cluster = dat$clust, R = dat$R_dd_wt)

all_res <- 
  data.frame(
    cluster = c(levels(clust)[clust],"$V^{CR2}$"),
    y = c(yi, NA),
    wt_full = res_wls_full,
    wt_absorb = res_wls_absorb,
    ols_full = res_ols_full,
    ols_absorb = res_ols_absorb,
    ols_full_het = res_ols_full_het,
    ols_absorb_het = res_ols_absorb_het
  )


# Validate against clubSandwich
library(clubSandwich)

# Ordinary least squares, ID working model

V_ols_full <- vcovCR(ols_fit, type = "CR2", cluster = dat$clust)
A_ols_full <- attr(V_ols_full,"adjustments")
emat_ols_full <- mapply(\(r, a) a %*% r, r = rddi_ols, a = A_ols_full)

V_ols_absorb <- vcovCR(ols_absorb, type = "CR2", cluster = dat$clust)
A_ols_absorb <- attr(V_ols_absorb, "adjustments")
emat_ols_absorb <- mapply(\(r, a) a %*% r, r = rddi_ols, a = A_ols_absorb)

# Ordinary least squares, heteroskedastic working model

V_ols_full_het <- vcovCR(ols_fit, type = "CR2", cluster = dat$clust, target = dat$R)
A_ols_full_het <- attr(V_ols_full_het,"adjustments")
emat_ols_full_het <- mapply(\(r, a) a %*% r, r = rddi_ols, a = A_ols_full_het)

V_ols_absorb_het <- vcovCR(ols_absorb, type = "CR2", cluster = dat$clust, target = dat$R)
A_ols_absorb_het <- attr(V_ols_absorb_het, "adjustments")
emat_ols_absorb_het <- mapply(\(r, a) a %*% r, r = rddi_ols, a = A_ols_absorb_het)

# Weighted least squares, heteroskedastic working model

V_wt_full <- vcovCR(wls_fit, type = "CR2", cluster = dat$clust, inverse_var = TRUE)
A_wt_full <- attr(V_wt_full, "adjustments")
emat_wt_full <- mapply(\(r, w, a) a %*% (w * r), r = rddi_wt, w = wi, a = A_wt_full)

V_wt_absorb <- vcovCR(wls_absorbed, type = "CR2", cluster = dat$clust, inverse_var = TRUE)
A_wt_absorb <- attr(V_wt_absorb, "adjustments")
emat_wt_absorb <- mapply(\(r, w, a) a %*% (w * r), r = rddi_wt, w = wi, a = A_wt_absorb)

RWA <- 
  data.frame(
    cluster = clust,
    y = yi,
    wt_full = unlist(emat_wt_full),
    wt_absorb = unlist(emat_wt_absorb),
    ols_full = unlist(emat_ols_full),
    ols_absorb = unlist(emat_ols_absorb),
    ols_full_het = unlist(emat_ols_full_het),
    ols_absorb_het = unlist(emat_ols_absorb_het)
  )

V_comp <- 
  data.frame(
    cluster = "$V^{CR2}$",
    y = NA,
    wt_full = V_wt_full["R","R"],
    wt_absorb = V_wt_absorb[1,1],
    ols_full = V_ols_full["R","R"],
    ols_absorb = V_ols_absorb[1,1],
    ols_full_het = V_ols_full_het["R","R"],
    ols_absorb_het = V_ols_absorb_het[1,1]
  )

all.equal(V_comp[,3:8], all_res[11,3:8], check.attributes = FALSE)
all.equal(RWA[,3:8], all_res[-11,3:8], check.attributes = FALSE)
