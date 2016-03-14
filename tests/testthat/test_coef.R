context("t-tests")

m <- 15
n <- 8 
cluster <- factor(rep(LETTERS[1:m], each = n))
N <- length(cluster)
X_btw <- rep(rep(LETTERS[1:3], c(3,5,7)), each = n)
X_wth <- rep(rep(c(0,1), each = n / 2), m)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
y <- nu + e

dat <- data.frame(y, X_btw, X_wth, cluster, row = 1:N)

