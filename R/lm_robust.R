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
#' ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)
#' 
#' if (requireNamespace("estimatr", quietly = TRUE)) withAutoprint({
#'   library(estimatr)
#' 
#'   lm_fit <- lm_robust(
#'     weight ~ Time + Diet:Time, 
#'     data = ChickWeight
#'    )
#'   vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#'   
#'   lm_fit_clust <- lm_robust(
#'     weight ~ Time + Diet:Time, data = ChickWeight,
#'     clusters = Chick
#'    )
#'   vcovCR(lm_fit_clust, type = "CR2")
#'   
#'   lm_fit_fe <- lm_robust(
#'     weight ~ Time:Diet, data = ChickWeight, 
#'     clusters = Chick, 
#'     fixed_effects = ~ Chick
#'    )
#'   vcovCR(lm_fit_fe, type = "CR2")
#' 
#'   if (requireNamespace("plm", quietly = TRUE)) withAutoprint({
#' 
#'     data("Produc", package = "plm")
#'     lm_individual <- lm_robust(
#'       log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#'       data = Produc,
#'       fixed_effects = ~ state,
#'       cluster = state
#'      )
#'     vcovCR(lm_individual, type = "CR2")
#'   
#'   })
#'   
#' })
#' 
#' @export

vcovCR.lm_robust <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  
  if (obj$fes && !requireNamespace("fixest", quietly = TRUE)) message("For improved performance in models with fixed effects, install the package {fixest}.")
  
  obj$model.frame <- model.frame(obj)
  
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
  
  model.frame(obj)[["(clusters)"]]
}


#' @export

augmented_model_matrix.lm_robust <- function(obj, cluster, inverse_var, ignore_FE) {
  
  if(!obj$fes) return(NULL)
  
  fe_formula <- as.formula(obj$call$fixed_effects)
  fe_formula <- update(fe_formula, ~ . - 1)
  mf <- model.frame(obj)
  model.matrix(fe_formula, data = mf)
  
}


requireNamespace <- function(...) base::requireNamespace(...)


#' @export
model_matrix.lm_robust <- function(obj) {
  
  if ("lm_lin" %in% as.character(obj$call[[1]])) {
    
    # Reconstruct the data and formulas
    data <- eval(obj$call$data, envir = environment(formula(obj)))
    
    frm <- as.formula(obj$call$formula)
    covariates <- as.formula(obj$call$covariates)
    
    if (is.null(covariates)) {
      # No covariates case - keep original formula as is
      # (it already has ~ 0 if that's what was specified)
      frm <- frm
    } else {
      # With covariates - preserve the intercept specification from original formula
      lhs <- deparse(frm[[2]])
      treatment <- all.vars(frm[[3]])
      covars <- all.vars(covariates)
      
      main_effects <- c(treatment, covars)
      interactions <- paste0(treatment, ":", covars)
      rhs_string <- paste(c(main_effects, interactions), collapse = " + ")
      
      # Check if original formula had no intercept (~ 0 or ~ -1)
      original_terms <- terms(frm)
      has_intercept <- attr(original_terms, "intercept") == 1
      
      if (has_intercept) {
        frm <- as.formula(paste(lhs, "~", rhs_string))
      } else {
        frm <- as.formula(paste(lhs, "~ 0 +", rhs_string))
      }
    }
    
    mf <- model.frame(frm, data)
    
    # Apply centering if needed
    if (!is.null(obj$scaled_center)) {
      covar_names <- names(obj$scaled_center)
      for (v in covar_names) {
        if (v %in% names(mf)) {
          mf[[v]] <- mf[[v]] - obj$scaled_center[[v]]
        }
      }
      
      Xmat <- model.matrix(frm, mf)
      
      # modify column names
      old_names <- colnames(Xmat)
      for (v in covar_names) {
        v_pattern <- paste0("(^",v,"$|:",v,")")
        i <- grepl(v_pattern, old_names)
        new_names <- gsub(v,paste0(v,"_c"), old_names[i])
        colnames(Xmat)[i] <- new_names
      }

    } else {
      Xmat <- model.matrix(frm, mf)
    }
    
    return(Xmat)
    
  }
  
  if (!obj$fes) {
    X_mat <- model.matrix(obj)
    
    return(X_mat)
  }
  
  # Fixed effects case (unchanged)
  mf <- model.frame(obj)
  
  frm <- as.formula(obj$call$formula)
  X_mat <- model.matrix(frm, data = mf)
  intercept_col <- colnames(X_mat) == "(Intercept)"
  if (any(intercept_col)) {
    X_mat <- X_mat[,!intercept_col,drop=FALSE]
  }
  
  fe_formula <- as.formula(obj$call$fixed_effects)
  
  if (requireNamespace("fixest", quietly = TRUE)) {
    fe_frame <- mf[attr(terms(fe_formula),"term.labels")]
    X_demean <- fixest::demean(X = X_mat, f = fe_frame)
  } else {
    fe_formula <- update(fe_formula, ~ . - 1)
    F_mat <- model.matrix(fe_formula, data = mf)
    X_reg <- stats::lm.fit(F_mat, X_mat)
    X_demean <- X_reg$residuals
  }
  
  X_demean
  
}


#' @export

model.frame.lm_robust <- function (formula, ...) {
  
  # Check if model.frame is already stored in the object
  mf <- formula$model.frame
  if (!is.null(mf)) return(mf)
  
  # Check if model was made using lm_lin
  if ("lm_lin" %in% as.character(formula$call[[1]])) {
    
    # Get original data
    data <- eval(formula$call$data, envir = environment(formula(formula)))
    
    # Extract the original formula and covariates from the call
    frm <- as.formula(formula$call$formula)
    covariates <- as.formula(formula$call$covariates)
    
    if (is.null(covariates)) {
      frm <- update(frm, ~ . - 1)
    } else {
      # Extract response and treatment
      lhs <- deparse(frm[[2]])
      treatment <- all.vars(frm[[3]])  # Assumes single treatment variable
      covars <- all.vars(covariates)
      
      # Main effects and interactions
      main_effects <- c(treatment, covars)
      interactions <- paste0(treatment, ":", covars)
      rhs_string <- paste(c(main_effects, interactions), collapse = " + ")
      
      # Construct full formula
      frm <- as.formula(paste(lhs, "~", rhs_string))
    }
    
    mf <- model.frame(frm, data)
    
    # Add (weights) column if weighted
    if (formula$weighted) {
      weights_expr <- formula$call$weights
      weights_val <- eval(weights_expr, envir = data, enclos = environment(formula(formula)))
      mf[["(weights)"]] <- weights_val
    }
    
    # Add (clusters) column if clustered
    if (formula$clustered) {
      cluster_expr <- formula$call$clusters
      cluster_val <- eval(cluster_expr, envir = data, enclos = environment(formula(formula)))
      mf[["(clusters)"]] <- cluster_val
    }
    
    # Center covariates if needed
    if (!is.null(formula$scaled_center)) {
      covar_names <- names(formula$scaled_center)
      for (v in covar_names) {
        if (v %in% names(mf)) {
          mf[[v]] <- mf[[v]] - formula$scaled_center[[v]]
        }
      }
    }
    
    return(mf)
  }
  
  
  # environment where initial call was evaluated
  fit_env <- environment(formula$terms) 
  
  # Extract relevant arguments from call
  cl <- formula$call
  mf_args <- match(c("formula","data","weights","subset","clusters"), names(cl), 0L)
  
  # Construct the model.frame for outcome and core predictors
  mf_cl <- cl[c(1L, mf_args)]
  mf_cl[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_cl, envir = fit_env)
  if (!formula$fes) return(mf)
  
  # Construct a model.frame for fixed effects
  fe_args <- match(c("fixed_effects","data","subset"), names(cl), 0L)
  fe_cl <- cl[c(1L, fe_args)]
  names(fe_cl)[[2]] <- "formula"
  fe_cl[[1L]] <- quote(stats::model.frame)
  mf_fe <- eval(fe_cl, envir = fit_env)
  
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
  
  mf <- model.frame(object)
  
  mf[[object$outcome]] - object$fitted.values
  
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

