



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
}

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

# written by GPT
model.frame.lm_robust <- function (object, ...) {
  # If a model frame is already stored, use it.
  if (!is.null(object$mf)) {
    return(object$mf)
  }
  # Otherwise, temporarily treat as an lm and extract the model frame.
  original_class <- class(object)
  class(object) <- "lm"
  mf <- model.frame(object, ...)
  # Optionally restore the class.
  class(object) <- original_class
  mf
}


residuals_CS.lm_robust <- function(obj) {
  # w <- obj$weights
  # if (is.null(w) || all(pos_wts <- w > 0)) {
  #   residuals(obj)
  # } else {
  #   residuals(obj)[pos_wts]
  # }
}