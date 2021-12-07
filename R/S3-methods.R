#----------------------------------------------
# get "working" variance-covariance matrix
#----------------------------------------------

targetVariance <- function(obj, cluster) UseMethod("targetVariance")

#' @export

targetVariance.default <- function(obj, cluster) {
  matrix_list(rep(1, length(cluster)), cluster, "both")
}

#----------------------------------------------
# get weighting matrix
#----------------------------------------------

weightMatrix <- function(obj, cluster) UseMethod("weightMatrix")

#' @export

weightMatrix.default <- function(obj, cluster) {
  weights <- weights(obj)
  if (is.null(weights)) {
    weights <- w_scale <- 1
  } else {
    w_scale <- mean(weights)
    weights <- weights / w_scale
  }
  W <- rep(weights, length.out = length(cluster))
  W_list <- matrix_list(W, cluster, "both")
  attr(W_list, "w_scale") <- w_scale
  W_list
}

#----------------------------------------------
# get X matrix
#----------------------------------------------

model_matrix <- function(obj) UseMethod("model_matrix")

#' @export

model_matrix.default <- function(obj) {
  model.matrix(obj)
}

#----------------------------------------------
# get augmented design matrix
#----------------------------------------------

augmented_model_matrix <- function(obj, cluster, inverse_var, ignore_FE) UseMethod("augmented_model_matrix")

#' @export

augmented_model_matrix.default <- function(obj, cluster, inverse_var, ignore_FE) {
  NULL
}

#----------------------------------------------
# get residuals
#----------------------------------------------

residuals_CS <- function(obj) UseMethod("residuals_CS") 

#' @export

residuals_CS.default <- function(obj) {
  residuals(obj)
}

#----------------------------------------------
# get coefficient estimates
#----------------------------------------------

coef_CS <- function(obj) UseMethod("coef_CS") 

#' @export

coef_CS.default <- function(obj) {
  coef(obj)
}

#----------------------------------------------
# get bread matrix
#----------------------------------------------

# bread matrices imported from sandwich package or elsewhere
#' @importFrom sandwich bread

get_bread <- function(obj) bread(obj)

v_scale <- function(obj) UseMethod("v_scale")

#' @export

v_scale.default <- function(obj) {
  nobs(obj)
}
