

# NOTE: has data field added, since it requires original data to get residuals
vcovCR.lm_robust <- function(obj, cluster, type, data, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) stop("You must specify a clustering variable.")
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  
  cluster <- droplevels(as.factor(cluster))
  
  X <- model_matrix(obj)
  if (any(alias)) {
    X <- X[, !alias, drop = FALSE]
  }  
  
  p <- NCOL(X)
  N <- NROW(X)
  
  cluster_length <- length(cluster)
  
  if (cluster_length != N) {
    
    cluster <- droplevels(handle_vectors(cluster, obj))
    
    if (length(cluster) != N) {
      stop("Clustering variable must have length equal to the number of rows in the data used to fit obj.")
    }
    
  } 
  
  if (any(is.na(cluster))) stop("Clustering variable cannot have missing values.")
  
  J <- nlevels(cluster)
  if (J < 2) stop("Cluster-robust variance estimation will not work when the data only includes a single cluster.")
  
  X_list <- matrix_list(X, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, w = W_list)
  
  if (is.null(target)) {
    if (inverse_var) {
      Theta_list <- lapply(W_list, function(w) chol2inv(chol(w)))
    } else {
      Theta_list <- targetVariance(obj, cluster)
    }
  } else {
    if (!is.list(target)) {
      if (length(target) != N) {
        target <- handle_vectors(target, obj)
      }
      Theta_list <- matrix_list(target, cluster, "both")
    } else {
      Theta_list <- target
    }
  }
  
  if (type %in% c("CR2","CR4")) {
    S <- augmented_model_matrix(obj, cluster, inverse_var, ignore_FE)
    
    if (is.null(S)) {
      rm(S)
      U_list <- X_list
      UW_list <- XW_list
    } else {
      U <- cbind(X, S)
      rm(S)
      U_list <- matrix_list(U, cluster, "row")
      UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
    }
    
    UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
    M_U <- matrix_power(Reduce("+",UWU_list), p = -1)
  }
  
  adjustments <- do.call(type, args = mget(names(formals(type))))
  
  E_list <- adjust_est_mats(type = type, est_mats = XW_list, adjustments = adjustments)
  
  resid <- residuals_CS(obj, data)
  res_list <- split(resid, cluster)
  
  components <- do.call(cbind, Map(function(e, r) e %*% r, e = E_list, r = res_list))
  
  v_scale <- v_scale(obj)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) w_scale <- 1L
  
  if (form == "estfun") {
    bread <- sandwich::bread(obj)
    estfun <- bread %*% components
    return(estfun * (w_scale / v_scale))
  }
  
  meat <- tcrossprod(components) * w_scale^2 / v_scale
  
  if (form == "sandwich") {
    bread <- sandwich::bread(obj)
  } else if (form == "meat") {
    bread <- NULL
  } else if (is.matrix(form)) {
    bread <- form
    form <- "sandwich"
  } 
  
  vcov <- switch(form, 
                 sandwich = bread %*% meat %*% bread / v_scale,
                 meat = meat)
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "type") <- type
  attr(vcov, "cluster") <- cluster
  attr(vcov, "bread") <- bread
  attr(vcov, "v_scale") <- v_scale
  attr(vcov, "est_mats") <- XW_list
  attr(vcov, "adjustments") <- adjustments
  attr(vcov, "target") <- Theta_list
  attr(vcov, "inverse_var") <- inverse_var
  attr(vcov, "ignore_FE") <- ignore_FE
  class(vcov) <- c("vcovCR","clubSandwich")
  return(vcov)
}


# model_matrix.lm_robust <- function(obj) {
  # model_matrix <- model.matrix(obj)
  # 
  # w <- obj$weights
  # if (is.null(w) || all(pos_wts <- w > 0)) {
  #   return(model_matrix)
  # } else {
  #   return(model_matrix[pos_wts > 0,,drop=FALSE])
# }

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
residuals_CS.lm_robust <- function(obj, data, ...) {
  w <- obj$weights
  if (is.null(w) || all(pos_wts <- w > 0)) {
    residuals(obj, data, ...)
  } else {
    residuals(obj, data, ...)[pos_wts]
  }
}

#' @export
residuals.lm_robust <- function(obj, data, ...) {
  col <- obj$terms[[2]]
  obj$fitted.values - data[[col]]
}