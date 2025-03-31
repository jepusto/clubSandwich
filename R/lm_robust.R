



vcovCR.lm_robust <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}


model_matrix.lm_robust <- function(obj) {
  # model_matrix <- model.matrix(obj)
  # 
  # w <- obj$weights
  # if (is.null(w) || all(pos_wts <- w > 0)) {
  #   return(model_matrix)
  # } else {
  #   return(model_matrix[pos_wts > 0,,drop=FALSE])
  # }
  
  
}