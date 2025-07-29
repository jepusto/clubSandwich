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
  
  # get formula
  frm <- as.formula(obj$call$formula)
  
  # If model was made using lm_lin
  if ("lm_lin" %in% as.character(obj$call[[1]])) {
    # get covariate formula, if it exists
    covariates <- as.formula(obj$call$covariates)
    
    if (!is.null(covariates)) { # No covariates case - keep original formula as is
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
    # get model frame
    mf <- model.frame(obj)
    # use model frame and formula built above to return model matrix
    return(model.matrix(frm, data = mf))
    
  } else { # If model was made using lm_robust
    # If no fixed effects, just return default mm
    if (!obj$fes) {
      return(model.matrix(obj))
    } else { # If fixed effects
      # get model frame
      mf <- model.frame(obj)
      
      # get base model matrix
      X_mat <- model.matrix(frm, data = mf)
      intercept_col <- colnames(X_mat) == "(Intercept)"
      if (any(intercept_col)) {
        X_mat <- X_mat[,!intercept_col,drop=FALSE]
      }
      
      # get fixed effects formula
      fe_formula <- as.formula(obj$call$fixed_effects)
      
      # use fixed effects formula to get model matrix for fixed effects,
      # then connect it to the base model matrix
      if (requireNamespace("fixest", quietly = TRUE)) {
        frame <- mf[attr(terms(fe_formula),"term.labels")]
        X_demean <- fixest::demean(X = X_mat, f = frame)
      } else {
        fe_formula <- update(fe_formula, ~ . - 1)
        F_mat <- model.matrix(fe_formula, data = mf)
        X_reg <- stats::lm.fit(F_mat, X_mat)
        X_demean <- X_reg$residuals
      }
      
      return(X_demean)
    }
  }
}


#' @export

model.frame.lm_robust <- function (formula, ...) {
  
  # Check if model.frame is already stored in the object
  mf <- formula$model.frame
  if (!is.null(mf)) return(mf)
  
  # environment where initial call was evaluated
  fit_env <- environment(formula$terms) 
  
  # Extract relevant arguments from call
  cl <- formula$call
  
  # Extract relevant arguments from call
  mf_args <- match(c("formula","data","weights","subset","clusters"), names(cl), 0L)
  
  # Construct the model.frame for outcome and core predictors
  mf_cl <- cl[c(1L, mf_args)]
  mf_cl[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_cl, envir = fit_env)
  
  # If model was made using lm_lin
  if ("lm_lin" %in% as.character(formula$call[[1]])) {
    word <- "covariates"
  } else { # If model was made using lm_robust
    if (!formula$fes) return(mf)
    word <- "fixed_effects"
  }
  
  # Construct a model.frame for covariates or fixed effects
  word_args <- match(c(word,"data","subset"), names(cl), 0L)
  word_cl <- cl[c(1L, word_args)]
  names(word_cl)[[2]] <- "formula"
  word_cl[[1L]] <- quote(stats::model.frame)
  mf_word <- eval(word_cl, envir = fit_env)
  # add centering here
  # compare omitted rows across model and covariates/fixed effects
  mf_omit <- na.action(mf)
  word_omit <- na.action(mf_word)
  
  # combine model.frames for model and for covariates/fixed effects  
  if (identical(mf_omit, word_omit)) {
    mf_combined <- cbind(mf, mf_word)
    mf_combined_omit <- mf_omit
  } else {
    mf_combined <- cbind(
      mf[!(rownames(mf) %in% word_omit),,drop=FALSE],
      mf_word[!(rownames(mf_word) %in% mf_omit),,drop=FALSE]
    )
    mf_combined_omit <- sort(c(mf_omit, word_omit))
    i_unique <- !duplicated(mf_combined_omit)
    mf_combined_omit <- mf_combined_omit[i_unique]
    class(mf_combined_omit) <- "omit"
  }
  
  # Apply centering if needed
  if (!is.null(formula$scaled_center)) {
    covar_names <- names(formula$scaled_center)
    for (v in covar_names) {
      if (v %in% names(mf_combined)) {
        mf_combined[[v]] <- mf_combined[[v]] - formula$scaled_center[[v]]
      }
    }
    
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

