#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an \code{estimatr::lm_robust} object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an 
#' \code{\link[estimatr]{lm_robust}} object.
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
#' data("ChickWeight", package = "datasets")
#' 
#' if (requireNamespace("estimatr", quietly = TRUE)) withAutoprint({
#' 
#'   lm_fit <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight)
#'   vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#'   
#'   lm_fit_clust <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight,
#'                             clusters = ChickWeight$Chick)
#'   vcovCR(lm_fit, type = "CR2")
#'   
#'   lm_fit_fe <- lm_robust(weight ~ 0 + Time:Diet, data = ChickWeight, 
#'   clusters = Chick, fixed_effects = ~Chick)
#'   vcovCR(lm_fit_fe, type = "CR2")
#' 
#'   if (requireNamespace("plm", quietly = TRUE)) withAutoprint({
#' 
#'     data("Produc", package = "plm")
#'     lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#'     individual_index <- !grepl("state", names(coef(lm_individual)))
#'     vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index]
#'   
#'   })
#'   
#' })
#' 
#' @export

vcovCR.lm_robust <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  
  if (missing(cluster)) {
    cluster <- findCluster.lm_robust(obj)
    if (is.null(cluster)) stop("You must specify a clustering variable or `obj` must include a clustering variable.")
  }
  
  if (is.null(inverse_var)) inverse_var <- is.null(weights(obj)) & is.null(target)
  
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}


#' Pulls clustering variable from lm_robust objects, if they have one.
#'
#' @param obj an lm_robust object
#'
#' @return The data within the clustering variable
#' @keywords internal
#' @noRd

findCluster.lm_robust <- function(obj) {
  
  if (!obj$clustered) return(NULL)
  
  cluster_expr <- obj$call$clusters
  fit_env <- environment(obj$terms)
  
  if (!is.null(obj$call$data)) {
    data_val <- eval(obj$call$data, envir = fit_env)
    cluster <- eval(cluster_expr, envir = data_val, enclos = fit_env)
  } else {
    cluster <- eval(cluster_expr, envir = fit_env)
  }
  
  return(cluster)
}


#' @export

augmented_model_matrix.lm_robust <- function(obj, cluster, inverse_var, ignore_FE) {
  
  if(!obj$fes) return(NULL)
  
  frm <- as.formula(obj$call$formula)
  
  fe_expr <- obj$call$fixed_effects[[2]]
  update_formula <- substitute(. ~ fe_expr + ., list(fe_expr = fe_expr))
  frm <- update(frm, update_formula)
  
  model.matrix(frm, model.frame(obj))
}


#' @export

model_matrix.lm_robust <- function(obj) {
  
  if (!obj$fes) return(model.matrix(obj))
  
  frm <- as.formula(obj$call$formula)
  X_mat <- model.matrix(frm, model.frame(obj))
  
  fe_frm <- as.formula(paste("~ 0 +", obj$call$fixed_effects[[2]]))
  
  F_mat <- model.matrix(fe_frm, model.frame(obj))
  
  model <- stats::lm.fit(F_mat, X_mat)
  
  model$residuals
  
}


#' @export

model.frame.lm_robust <- function (formula, ...) {
  
  # Extract relevant arguments
  cl <- formula$call
  m <- match(c("formula","data","weights","subset","clusters"), names(cl), 0L)
  mf <- cl[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)

  # check for fixed_effects
  if (formula$fes) {
    fe_expr <- formula$call$fixed_effects[[2]]
    update_formula <- substitute(. ~ fe_expr + ., list(fe_exp = fe_expr))
    mf$formula <- update(eval(mf$formula), update_formula)
  }

  # Construct the model.frame
  mf <- eval(mf, parent.frame())

  mf
}



#' @export

residuals.lm_robust <- function(object, ...) {
  
  model.frame(object)[[object$outcome]] - object$fitted.values
  
}


#' @export

bread.lm_robust <- function(x, ...) {
  
  N <- nobs(x)
  
  X_mat <- model_matrix(x)
  
  if (x$weighted) {
    XtWX <- crossprod(X_mat, x$weights * X_mat)
  } else {
    XtWX <- crossprod(X_mat)
  }
  
  N * solve(XtWX)
}


#' @export

na.action.lm_robust <- function(object, ...)  {
  
  na.action(model.frame(object))
  
}

