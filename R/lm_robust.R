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

compare_model_frames <- function(data1, data2, ...) {
  cl <- match.call()
  
  cl1 <- cl
  cl1$data2 <- NULL
  names(cl1)[2] <- "data"
  cl1[[1]] <- quote(estimatr::lm_robust)
  m1 <- suppressWarnings(eval(cl1, parent.frame()))
  mf1 <- model.frame(m1)

  cl2 <- cl
  cl2$data1 <- NULL
  names(cl2)[2] <- "data"
  cl2[[1]] <- quote(estimatr::lm_robust)
  m2 <- suppressWarnings(eval(cl2, parent.frame()))
  mf2 <- model.frame(m2)
  
  testthat::expect_equivalent(mf1, mf2)
}

#' @export

model.frame.lm_robust <- function (formula, ...) {
  
  # Extract relevant arguments from call
  cl <- formula$call
  mf_args <- match(c("formula","data","weights","subset","clusters"), names(cl), 0L)
  
  # Construct the model.frame
  mf_cl <- cl[c(1L, mf_args)]
  mf_cl[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_cl, parent.frame())
  if (!formula$fes) return(mf)
  
  # Construct a model.frame for fixed effects
  fe_args <- match(c("fixed_effects","data","subset"), names(cl), 0L)
  fe_cl <- cl[c(1L, fe_args)]
  names(fe_cl)[[2]] <- "formula"
  fe_cl[[1L]] <- quote(stats::model.frame)
  mf_fe <- eval(fe_cl, parent.frame())
  
  # compare omitted rows across model and fixed effects
  mf_omit <- na.action(mf)
  fe_omit <- na.action(mf_fe)
  
  # combine model.frames for model and for fixed effects  
  if (identical(mf_omit, fe_omit)) {
    mf_combined <- cbind(mf, mf_fe)
    mf_combined_omit <- mf_omit
  } else {
    mf_combined <- cbind(
      mf[!(rownames(mf) %in% fe_omit),,drop=FALSE],
      mf_fe[!(rownames(mf_fe) %in% mf_omit),,drop=FALSE]
    )
    mf_combined_omit <- sort(c(mf_omit, fe_omit))
    i_unique <- !duplicated(mf_combined_omit)
    mf_combined_omit <- mf_combined_omit[i_unique]
    class(mf_combined_omit) <- "omit"
  }
  attr(mf_combined, "terms") <- attr(mf, "terms")
  attr(mf_combined, "na.action") <- mf_combined_omit
  
  return(mf_combined)
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

