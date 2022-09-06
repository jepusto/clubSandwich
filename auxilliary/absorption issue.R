library(lme4)
set.seed(20220906)
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
bread_inv <- 
  Map(function(x, w) t(x) %*% w %*% x, x = U_dot, w = W_list) %>%
  Reduce("+", .) 

bread <- chol2inv(chol(bread_inv))

beta <- 
  Map(function(x, w, y) t(x) %*% w %*% y, x = U_dot, w = W_list, y = y_dot) %>%
  Reduce("+", .) %>%
  solve(bread_inv, .)

as.vector(beta)
fixef(lme_dummy)[c("U.1","U.2")]

lme_absorb <- lmer(y ~ 0 + U.1 + U.2 + (0 + U.1 + U.2 || id), 
                   data = dat_absorb,
                   start = list(getME(lme_dummy, "theta")))
fixef(lme_absorb)
A_dummy <- attr(vcovCR(lme_dummy, type = "CR2"), "adjustments")
A_absorb <- attr(vcovCR(lme_absorb, type = "CR2"), "adjustments")

UtWA_dummy <- Map(function(u, w, a) t(u) %*% w %*% a, 
                  u = U_dot, w = W_list, a = A_dummy)
UtWA_absorb <- Map(function(u, w, a) t(u) %*% w %*% a, 
                   u = U_dot, w = W_list, a = A_absorb)

UtWA_dummy[[1]]
UtWA_absorb[[1]]
