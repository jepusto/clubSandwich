#----------------------------------------------
# get "working" variance-covariance matrix
#----------------------------------------------

targetVariance <- function(obj) UseMethod("targetVariance")

#----------------------------------------------
# get weighting matrix
#----------------------------------------------

weightMatrix <- function(obj) UseMethod("weightMatrix")

#----------------------------------------------
# get X matrix
#----------------------------------------------

model_matrix <- function(obj) UseMethod("model_matrix")

model_matrix.default <- function(obj) {
  model.matrix(obj)
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

