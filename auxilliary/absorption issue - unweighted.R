set.seed(20220926)
m <- 4                                             # number of clusters
ni <- 2 + rpois(m, 3.5)                            # cluster sizes
N <- sum(ni)                                       # total sample size
id <- factor(rep(LETTERS[1:m], ni))                # cluster ID
R <- rnorm(N)                                      # focal predictor
dat <- data.frame(R, id)                           # create raw data frame
X <- model.matrix(~ R + id + 0, data = dat)        # full predictor matrix
Ui <- tapply(R, id, \(x) x - mean(x))              # absorbed version of R
U <- unsplit(Ui, id)
matrix_power <- function(x, p) {
  eig <- eigen(x, symmetric = TRUE)
  val_p <- with(eig, ifelse(values > 10^-12, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

MX <- solve(crossprod(X))
B <- 
  by(X, id, as.matrix) |>
  lapply(\(x) diag(nrow(x)) - x %*% MX %*% t(x))
A <- lapply(B, matrix_power, p = -1/2)

MU <- 1 / crossprod(U)
Btilde <- lapply(Ui, \(x) diag(length(x)) - x %*% MU %*% t(x))
Atilde <- lapply(Btilde, matrix_power, p = -1/2)

all.equal(A, Atilde)

UiAtilde <- mapply(\(u, a) t(u) %*% a, u = Ui, a = Atilde, SIMPLIFY = FALSE)
UiA <- mapply(\(u, a) t(u) %*% a, u = Ui, a = A, SIMPLIFY = FALSE)
all.equal(UiAtilde, UiA)

Tmat <- model.matrix(~ 0 + id)
Ti <- by(Tmat, id, as.matrix)
MT <- solve(crossprod(Tmat))
tiMTtit <- lapply(Ti, \(x) x %*% MT %*% t(x))
B_alt <- mapply(`-`, Btilde, tiMTtit)
all.equal(B, B_alt)

mapply(\(a, t) a %*% t, a = Atilde, t = Ti)
