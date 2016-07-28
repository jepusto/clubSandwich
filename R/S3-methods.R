#----------------------------------------------
# get "working" variance-covariance matrix
#----------------------------------------------

targetVariance <- function(obj, cluster) UseMethod("targetVariance")

targetVariance.default <- function(obj, cluster) {
  matrix_list(rep(1, nobs(obj)), cluster, "both")
}

#----------------------------------------------
# get weighting matrix
#----------------------------------------------

weightMatrix <- function(obj, cluster) UseMethod("weightMatrix")

weightMatrix.default <- function(obj, cluster) {
  weights <- weights(obj)
  if (is.null(weights)) weights <- 1
  W <- rep(weights, length.out = nobs(obj))
  matrix_list(W, cluster, "both")
}

#----------------------------------------------
# get X matrix
#----------------------------------------------

model_matrix <- function(obj) UseMethod("model_matrix")

model_matrix.default <- function(obj) {
  model.matrix(obj)
}

#----------------------------------------------
# get projection matrix
#----------------------------------------------

projection_matrix <- function(obj) UseMethod("projection_matrix")

projection_matrix.default <- function(obj) {
  model_matrix(obj)
}

#----------------------------------------------
# get augmented design matrix
#----------------------------------------------

augmented_model_matrix <- function(obj, cluster, inverse_var) UseMethod("augmented_model_matrix")

augmented_model_matrix.default <- function(obj, cluster, inverse_var) {
  NULL
}

#----------------------------------------------
# get residuals
#----------------------------------------------

residuals_CS <- function(obj) UseMethod("residuals_CS") 

residuals_CS.default <- function(obj) {
  residuals(obj)
}

#----------------------------------------------
# get coefficient estimates
#----------------------------------------------

coef_CS <- function(obj) UseMethod("coef_CS") 

coef_CS.default <- function(obj) {
  coef(obj)
}

