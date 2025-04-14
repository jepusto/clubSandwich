

#' @export
# vcovCR.lm_robust <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  
# }


# model_matrix.lm_robust <- function(obj) {
  # model_matrix <- model.matrix(obj)
  #
  # w <- obj$weights
  # if (is.null(w) || all(pos_wts <- w > 0)) {
  #   return(model_matrix)
  # } else {
  #   return(model_matrix[pos_wts > 0,,drop=FALSE])
# }
#' @export
model.matrix.lm_robust <- function (object, ...) {
  if (n_match <- match("x", names(object), 0L))
    object[[n_match]]
  else {
    data <- model.frame(object, ...)
    if (exists(".GenericCallEnv", inherits = FALSE))
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call("model.matrix.default", c(list(object = object,
                                             data = data, contrasts.arg = object$contrasts),
                                        dots))
    }
  }
}

# written by GPT, slightly edited by Sam
#' @export
model.frame.lm_robust <- function (obj, ...) {
  # If a model frame is already stored, use it.
  if (!is.null(obj$mf)) {
    return(obj$mf)
  }
  # Otherwise, temporarily treat as an lm and extract the model frame.
  original_class <- class(obj)
  class(obj) <- "lm"
  mf <- model.frame(obj, ...)
  # Optionally restore the class.
  class(obj) <- original_class
  mf
}

#' @export
residuals.lm_robust <- function(obj, ...) {
  # data <- eval(obj$call$data, envir = parent.frame())
  # col <- obj$outcome
  # data[[col]] - obj$fitted.values
  model.frame(obj)[[obj$outcome]] - obj$fitted.values # from github discussion
}


#' @export
bread.lm_robust <- function(obj, ...) {
  
  N <- nobs(obj)
  
  X <- model_matrix(obj)
  
  w <- model.frame(obj)$weight
  
  XtWX <- crossprod(X, w * X)
  
  # p <- obj$rank
  # p1 <- 1L:p
  # R <- chol2inv(obj$qr[p1, p1, drop = FALSE])
  # return(N * R)
  
  return(N * solve(XtWX))
}




