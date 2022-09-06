
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

# turn block-diagonal into regular matrix

unblock <- function(A, block = attr(A, "groups")) {
  
  if (is.null(block)) block <- factor(rep(names(A), times = sapply(A, function(x) dim(x)[1])))
  n <- length(block)
  mat <- matrix(0, n, n)
  for (i in levels(block)) {
    index <- i == block
    mat[index,index] <- A[[i]]
  }
  return(mat)
}

matrix_power <- function(x, p, symmetric = TRUE, tol = -12) {
  eig <- eigen(x, symmetric = symmetric)
  val_p <- with(eig, ifelse(values > 10^tol, values^p, 0))
  with(eig, vectors %*% (val_p * t(vectors)))
}

chol_psd <- function(x) with(eigen(x, symmetric=TRUE), sqrt(pmax(values,0)) * t(vectors))


add_submatrices <- function(indices, small_mat, big_mat) {
  levs <- levels(indices)
  for (i in 1:length(levs)) {
    ind <- levs[i] == indices
    big_mat[ind,ind] <- small_mat[[i]] + big_mat[ind,ind]
  }
  big_mat
}
  
add_bdiag <- function(small_mats, big_mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  small_mats <- split(small_mats, big_indices)
  Map(add_submatrices, indices = small_indices, small_mat = small_mats, big_mat = big_mats)
}

nest_bdiag <- function(mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  mat_groups <- split(mats, big_indices)
  Map(unblock, A = mat_groups)
}

