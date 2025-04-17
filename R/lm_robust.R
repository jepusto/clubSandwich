#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an estimatr::lm_robust object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an 
#' \code{\link{estimatr::lm_robust}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{estimatr::lm_robust} objects.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be an identity matrix.
#' @inheritParams vcovCR
#'   
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists
#'   of a matrix of the estimated variance of and covariances between the
#'   regression coefficient estimates.
#'   
#' @seealso \code{\link{vcovCR}}
#' 
#' @examples 
#' 
#' NOTE: These are (currently) the same as the examples from vcovCR.lm
#' 
#' data("ChickWeight", package = "datasets")
#' lm_fit <- lm(weight ~ Time + Diet:Time, data = ChickWeight)
#' vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#' 
#' if (requireNamespace("plm", quietly = TRUE)) withAutoprint({
#' 
#'   data("Produc", package = "plm")
#'   lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#'   individual_index <- !grepl("state", names(coef(lm_individual)))
#'   vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index]
#' 
#'   # compare to plm()
#'   plm_FE <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#'                      data = Produc, index = c("state","year"), 
#'                      effect = "individual", model = "within")
#'   vcovCR(plm_FE, type="CR2")
#'   
#' })
#' 
#' @export
vcovCR.lm_robust <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) {
    if (obj$clustered) cluster <- get_cluster(obj) # pull clustering var from obj
    else stop("You must specify a clustering variable or your object must have one.")
  }
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}


#' Helper function written by GPT, edited by Sam
get_cluster <- function(obj) {
  # 1. Grab the unevaluated clusters argument
  cluster_expr <- obj$call$clusters
  
  # 2. Use the formula/environment of the fit as our evaluation backbone
  fit_env <- environment(obj$terms)
  
  # 3. If the user passed a data= argument, pull that data in...
  if (!is.null(obj$call$data)) {
    data_val <- eval(obj$call$data, envir = fit_env)
    # ...and evaluate the clusters expression inside that data
    cluster <- eval(cluster_expr, envir = data_val, enclos = fit_env)
  } else {
    # otherwise just eval it in the fit’s environment
    cluster <- eval(cluster_expr, envir = fit_env)
  }
  
  # 4. Done
  cluster
}


#' Same as model.matrix.lm
#' @export
model.matrix.lm_robust <- function (object, ...) 
{
  if (n_match <- match("x", names(object), 0L)) 
    object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
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
  
  if(obj$weighted) {
    XtWX <- crossprod(X, obj$weights * X)
  }
  else {
    XtWX <- crossprod(X)
  }
  
  return(N * solve(XtWX))
}




