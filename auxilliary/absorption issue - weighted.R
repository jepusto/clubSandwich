library(lme4)
set.seed(20220914)
J <- 5
nj <- 2 + rpois(J, 7)
N <- sum(nj)

id <- factor(rep(LETTERS[1:J], nj))
tm <- unsplit(lapply(nj, \(x) 1:x), id)
bU <- rnorm(J, mean = 1)[id] + rnorm(J, mean = 0.25, sd = 0.25)[id] * tm
Umat <- matrix(rnorm(2 * N, mean = bU), N, 2)
U_list <- matrix_list(Umat, id, "row")
T_list <- matrix_list(cbind(1, tm), id, "row")
T_full <- matrix_list(model.matrix(~ 0 + id + id:tm), id, "row")

bT <- rnorm(J, mean = 1)[id] + rnorm(J, mean = 0.25, sd = 0.25)[id] * tm
y <- rowSums(Umat) + bT + rnorm(N)
y_list <- matrix_list(as.matrix(y), id, "row")

dat <- cbind(data.frame(y, id, tm), U = Umat)

lme_dummy <- lmer(y ~ 0 + id + id:tm + U.1 + U.2 + (0 + U.1 + U.2 || id), data = dat)
W_list <- weightMatrix(lme_dummy)
Phi_list <- targetVariance(lme_dummy, cluster = id)
D_list <- lapply(Phi_list, chol)
D_inv <- lapply(D_list, solve)

# Absorbed matrices
reg_by <- function(x, w, y) {
  lhs <- t(x) %*% w %*% x
  rhs <- t(x) %*% w %*% y
  beta <- solve(lhs, rhs)
  y - x %*% beta
}
U_dot <- mapply(reg_by, x = T_list, w = W_list, y = U_list)
y_dot <- mapply(reg_by, x = T_list, w = W_list, y = y_list)
dat_absorb <- data.frame(U = do.call(rbind, U_dot), y = do.call(rbind, y_dot), id)

# Should be 2X2 matrices of zeros
mapply(function(x, w, y) t(x) %*% w %*% y, 
       x = T_list, w = W_list, y = U_dot, SIMPLIFY = FALSE)

# Calculate absorbed beta coefficients
UtWU <- 
  Map(function(x, w) t(x) %*% w %*% x, x = U_dot, w = W_list) %>%
  Reduce("+", .) 

MUd <- chol2inv(chol(UtWU))

beta <- 
  Map(function(x, w, y) t(x) %*% w %*% y, x = U_dot, w = W_list, y = y_dot) %>%
  Reduce("+", .) %>%
  solve(UtWU, .)

as.vector(beta)
fixef(lme_dummy)[c("U.1","U.2")]

# Calculate adjustment matrices by hand

TtWT_list <- Map(function(x, w) t(x) %*% w %*% x, x = T_full, w = W_list)
TtWT <- Reduce("+", TtWT_list) 
MT <- chol2inv(chol(TtWT))
TMTTt <- Map(\(x) x %*% MT %*% t(x), x = T_full)

B_list <- Map(function(d, phi, u, tmat) d %*% (phi - u %*% MUd %*% t(u) - tmat %*% MT %*% t(tmat)) %*% t(d), 
                  d = D_list, phi = Phi_list, u = U_dot, tmat = T_full)
B_ginv <- Map(matrix_power, B_list, p = -1)

Btilde_list <- Map(function(d, phi, u, tmat) d %*% (phi - u %*% MUd %*% t(u)) %*% t(d), 
                    d = D_list, phi = Phi_list, u = U_dot, tmat = T_full)
Btilde_ginv <- Map(matrix_power, Btilde_list, p = -1)

B_diff <- Map(\(x,y) x - y, x = Btilde_ginv, y = B_ginv)
B_diff

A_list <- Map(\(b, d) t(d) %*% matrix_power(b, p = -1/2) %*% d, b = B_list, d = D_list)
Atilde_list <- Map(\(b, d) t(d) %*% matrix_power(b, p = -1/2) %*% d, b = Btilde_list, d = D_list)
UtWA <- Map(\(u,w,a) t(u) %*% w %*% a, u = U_dot, w = W_list, a = A_list)
UtWAtilde <- Map(\(u,w,a) t(u) %*% w %*% a, u = U_dot, w = W_list, a = Atilde_list)
UtWA[[4]]
UtWAtilde[[4]]

Psi_list <- Map(function(phi, u) matrix_power(phi - u %*% MUd %*% t(u), -1), 
                phi = Phi_list, u = U_dot)
Psi_T <- Map(function(psi, tmat) psi %*% tmat, psi = Psi_list, tmat = T_list)
W_T <- Map(function(w, tmat) w %*% tmat, w = W_list, tmat = T_list)
all.equal(Psi_T, W_T)

MTspace_list <- lapply(TtWT_list, \(twt) MT - MT %*% twt %*% MT)
MTspace_inv <- lapply(MTspace_list, matrix_power, p = -1)
MTspace_check <- Map(function(a, ai) a %*% ai %*% a, a = MTspace_list, ai = MTspace_inv)
all.equal(MTspace_list, MTspace_check)

Tspace_list <- Map(function(w, tmat, mtspace) w %*% tmat %*% MT %*% mtspace %*% MT %*% t(tmat) %*% w,
                   w = W_list, tmat = T_full, mtspace = MTspace_list)
  
A_dummy <- attr(vcovCR(lme_dummy, type = "CR2"), "adjustments")
all.equal(A_dummy[[1]], A_list[[1]], check.attributes = FALSE)
all.equal(A_dummy[[2]], A_list[[2]], check.attributes = FALSE)
all.equal(A_dummy[[3]], A_list[[3]], check.attributes = FALSE)
all.equal(A_dummy[[4]], A_list[[4]], check.attributes = FALSE)
all.equal(A_dummy[[5]], A_list[[5]], check.attributes = FALSE)
