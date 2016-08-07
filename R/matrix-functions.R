
#---------------------------------------------
# matrix manipulation functions
#---------------------------------------------

sub_f <- function(x, fac, dim) {
  function(f) switch(dim,
                     row = x[fac==f, ,drop=FALSE],
                     col = x[ ,fac==f, drop=FALSE],
                     both = x[fac==f, fac==f, drop=FALSE])
}

matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") stop(paste0("Object must be a matrix in order to subset by ",dim,"."))
    x_list <- split(x, fac)
    lapply(x_list, function(x) diag(x, nrow = length(x)))
  } else {
    lapply(levels(fac), sub_f(x, fac, dim)) 
  }
}

matrix_power <- function(x, p, symmetric = TRUE, tol = -12) {
  eig <- eigen(x, symmetric = symmetric)
  val_p <- with(eig, ifelse(values > 10^tol, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

chol_psd <- function(x) with(eigen(x, symmetric=TRUE), sqrt(pmax(values,0)) * t(vectors))

