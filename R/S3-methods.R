#----------------------------------------------
# get "working" variance-covariance matrix
#----------------------------------------------

targetVariance <- function(obj, cluster) UseMethod("targetVariance")

#----------------------------------------------
# get weighting matrix
#----------------------------------------------

weightMatrix <- function(obj, cluster) UseMethod("weightMatrix")

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

residuals_CR <- function(obj) UseMethod("residuals_CR") 

residuals_CR.default <- function(obj) {
  residuals(obj)
}

#----------------------------------------------
# get coefficient estimates
#----------------------------------------------

coef_CR <- function(obj) UseMethod("coef_CR") 

coef_CR.default <- function(obj) {
  coef(obj)
}

