#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for an \code{estimatr::lm_robust}
#' object.
#'
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
#' of a set of regression coefficient estimates from an
#' \code{\link[estimatr]{lm_robust}} object.
#'
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. If not specified, will be detected from the
#'   \code{clusters} argument of \code{obj}.
#' @param type Character string specifying which small-sample adjustment should
#'   be used, with available options \code{"CR0"}, \code{"CR1"}, \code{"CR1p"},
#'   \code{"CR1S"}, \code{"CR2"}, or \code{"CR3"}. If not specified, will be
#'   detected from the \code{se_type} argument of \code{obj}. See "Details"
#'   section of \code{\link{vcovCR}} for further information.
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
#'   conf_int(lm_fit_clust, vcov = "CR2")
#'
#'   # similar model via lm_lin()
#'   lin_fit_clust <- lm_lin(
#'     weight ~ Diet, 
#'     covariates = ~ Time,
#'     data = ChickWeight,
#'     clusters = Chick
#'   )
#'   conf_int(lin_fit_clust, vcov = "CR2")
#'   
#'   lm_fit_fe <- lm_robust(
#'     weight ~ Time:Diet, data = ChickWeight,
#'     clusters = Chick,
#'     fixed_effects = ~ Chick
#'    )
#'   vcovCR(lm_fit_fe)
#'   
#'   # two-way fixed effects model
#'   data("MortalityRates")
#'   MortalityRates <- subset(MortalityRates, cause == "Motor Vehicle")
#'   MortalityRates$state <- factor(MortalityRates$state)
#'   MortalityRates$year <- factor(MortalityRates$year)
#'   MLDA_fit <- lm_robust(
#'     mrate ~ legal + beertaxa + beerpercap + winepercap + spiritpercap,
#'     fixed_effects = ~ year + state,
#'     data = MortalityRates,
#'     cluster = state
#'   )
#'   conf_int(MLDA_fit, vcov = "CR2")
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
  
  if (obj$fes && !fixest_is_available()) {
    message("For improved performance in models with fixed effects, install the package {fixest}.")
  }
  
  obj$model.frame <- model.frame(obj)
  
  if (missing(cluster)) {
    cluster <- findCluster.lm_robust(obj)
    if (is.null(cluster)) stop("You must specify a clustering variable or `obj` must include a clustering variable.")
  }
  
  if (missing(type)) {
    type <- switch(obj$se_type, CR0 = "CR0", CR2 = "CR2", stata = "CR1S", "No valid SE type")
    if (type == "No valid SE type") stop("You must specify a `type` of sandwich estimator to calculate or `obj` must include an `se_type` of 'CR0','CR2',or 'stata'.")
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
  
  model.frame(obj)[["(cluster)"]]
}


#' @export

augmented_model_matrix.lm_robust <- function(obj, cluster, inverse_var, ignore_FE) {
  
  if(!obj$fes) return(NULL)
  
  # get formula for the fixed effects
  fe_formula <- as.formula(obj$call$fixed_effects)
  fe_formula <- update(fe_formula, ~ . - 1)
  
  # get the model.frame
  mf <- model.frame(obj)
  
  # ensure fixed effects are all factors
  varnames <- all.vars(fe_formula)
  for (v in varnames) mf[[v]] <- as.factor(mf[[v]])
  
  # compute model.matrix of the fixed effects
  model.matrix(fe_formula, data = mf)
  
}


#' @export

model_matrix.lm_robust <- function(obj) {
  
  # get formula
  frm <- as.formula(obj$call$formula)
  
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
    if (fixest_is_available()) {
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


fixest_is_available <- function() requireNamespace("fixest", quietly = TRUE)


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

