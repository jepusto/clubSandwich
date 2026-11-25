#' Cluster-robust variance-covariance matrix for a glmmTMB object
#'
#' Calculates a cluster-robust variance estimator for the conditional
#' fixed effects of a Gaussian identity-link model fitted by
#' \code{\link[glmmTMB]{glmmTMB}}.
#'
#' @param cluster Optional vector identifying independent clusters. If omitted,
#'   an outermost random-effect grouping factor is used, provided that all
#'   other random effects are nested within it. Otherwise it must be specified.
#' @param target Optional working variance-covariance model for the small-sample
#'   adjustments. By default, the fitted marginal covariance is used.
#' @inheritParams vcovCR
#'
#' @details
#' The method supports Gaussian identity-link models without zero inflation,
#' prior weights, priors, dispersion random effects, or mapped conditional
#' fixed effects. Random-effect covariance structures are currently limited to
#' unstructured and diagonal structures. Every random-effect grouping factor
#' must be nested within the supplied clusters. Fixed-effects dispersion
#' models and offsets are supported.
#'
#' @return An object of class \code{c("vcovCR", "clubSandwich")}.
#' @seealso \code{\link{vcovCR}}
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE) &&
#'     requireNamespace("lme4", quietly = TRUE)) {
#'   data("sleepstudy", package = "lme4")
#'   fit <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
#'   vcovCR(fit, type = "CR2")
#'   coef_test(fit, vcov = "CR2")
#' }
#' @export
vcovCR.glmmTMB <- function(obj, cluster, type, target, inverse_var, form = "sandwich", ...) {
  validate_glmmTMB(obj)
  groups <- glmmTMB_groups(obj)
  if (missing(cluster)) {
    if (!length(groups)) {
      stop("A clustering variable must be supplied for models without random effects.")
    }
    candidates <- order(vapply(groups, function(g) nlevels(droplevels(g)), integer(1)))
    outer <- candidates[vapply(candidates, function(i) {
      all(vapply(groups, check_nested, outer_grp = groups[[i]], FUN.VALUE = logical(1)))
    }, logical(1))]
    if (!length(outer)) {
      stop("No outermost random-effect grouping factor; supply a clustering variable.")
    }
    cluster <- groups[[outer[1]]]
  }

  if (length(cluster) != nrow(model.frame(obj))) {
    omitted <- attr(model.frame(obj), "na.action")
    if (inherits(omitted, c("omit", "exclude"))) cluster <- cluster[-omitted]
  }
  if (length(cluster) != nrow(model.frame(obj))) {
    stop("Clustering variable must have length equal to the number of rows in the data used to fit obj.")
  }
  if (anyNA(cluster)) stop("Clustering variable cannot have missing values.")
  if (!all(vapply(groups, check_nested, outer_grp = cluster, FUN.VALUE = logical(1)))) {
    stop("Random effects must be nested within the clustering variable.")
  }

  if (missing(target)) target <- NULL
  if (missing(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, target = target,
          inverse_var = inverse_var, form = form)
}

validate_glmmTMB <- function(obj) {
  fam <- family(obj)
  if (fam$family != "gaussian" || fam$link != "identity") {
    stop("Only Gaussian models with an identity link are supported.")
  }
  if (!is.null(weights(obj))) stop("Models with prior weights are not currently supported.")
  if (length(glmmTMB::fixef(obj)$zi) || ncol(obj$obj$env$data$Zzi) > 0) {
    stop("Zero-inflated models are not currently supported.")
  }
  if (length(obj$modelInfo$reStruc$dispReStruc)) {
    stop("Models with random effects in the dispersion model are not currently supported.")
  }
  if (!is.null(obj$modelInfo$map$beta)) {
    stop("Models with mapped conditional fixed effects are not currently supported.")
  }
  if (!is.null(obj$modelInfo$priors) && nrow(obj$modelInfo$priors) > 0) {
    stop("Models fitted with priors are not currently supported.")
  }
  terms <- obj$modelInfo$reStruc$condReStruc
  codes <- vapply(terms, function(x) names(x$blockCode), character(1))
  if (any(!codes %in% c("us", "diag", "homdiag"))) {
    stop("Only unstructured and diagonal random-effect covariance structures are currently supported.")
  }
}

glmmTMB_groups <- function(obj) {
  lapply(obj$modelInfo$reTrms$cond$flist, function(g) droplevels(as.factor(g)))
}

# Connected components of random-effect groupings give independent covariance blocks.
glmmTMB_components <- function(obj) {
  groups <- glmmTMB_groups(obj)
  component <- seq_len(nrow(model.frame(obj)))
  repeat {
    previous <- component
    for (g in groups) {
      group_min <- tapply(component, g, min)
      component <- pmin(component, unname(group_min[as.character(g)]))
    }
    if (identical(component, previous)) break
  }
  factor(component)
}

#' @export
coef_CS.glmmTMB <- function(obj) glmmTMB::fixef(obj)$cond

#' @export
model_matrix.glmmTMB <- function(obj) {
  as.matrix(model.matrix(obj, component = "cond", include_rankdef = TRUE))
}

#' @export
residuals_CS.glmmTMB <- function(obj) {
  beta <- coef_CS(obj)
  X <- model_matrix(obj)[, !is.na(beta), drop = FALSE]
  y <- model.response(model.frame(obj))
  offset <- obj$obj$env$data$offset
  as.numeric(y - X %*% beta[!is.na(beta)] - offset)
}

#' @export
targetVariance.glmmTMB <- function(obj, cluster = glmmTMB_components(obj)) {
  validate_glmmTMB(obj)
  cluster <- droplevels(as.factor(cluster))
  terms <- obj$modelInfo$reStruc$condReStruc
  sd <- as.numeric(predict(obj, type = "disp"))
  omitted <- attr(model.frame(obj), "na.action")
  if (length(sd) != length(cluster) && inherits(omitted, "exclude") &&
      length(sd) == length(cluster) + length(omitted)) {
    sd <- sd[-omitted]
  }
  if (length(sd) != length(cluster) || any(!is.finite(sd)) || any(sd <= 0)) {
    stop("The fitted residual standard deviations must be finite and positive.")
  }

  Z <- lme4::getME(obj, "Z")
  if (length(terms)) {
    covariances <- glmmTMB::VarCorr(obj)$cond
    if (length(covariances) != length(terms) ||
        any(vapply(seq_along(terms), function(i) {
          !identical(dim(covariances[[i]]), rep(as.integer(terms[[i]]$blockSize), 2))
        }, logical(1)))) {
      stop("Could not reconstruct the random-effect covariance matrices.")
    }
    blocks <- unlist(lapply(seq_along(terms), function(i) {
      rep(list(as.matrix(covariances[[i]])), terms[[i]]$blockReps)
    }), recursive = FALSE)
    G <- do.call(Matrix::bdiag, blocks)
  }

  rows <- split(seq_along(cluster), cluster)
  lapply(rows, function(i) {
    V <- diag(sd[i]^2, nrow = length(i))
    if (length(terms)) {
      Zi <- Z[i, , drop = FALSE]
      V <- V + as.matrix(Zi %*% G %*% Matrix::t(Zi))
    }
    V
  })
}

#' @export
weightMatrix.glmmTMB <- function(obj, cluster = glmmTMB_components(obj)) {
  lapply(targetVariance(obj, cluster), function(v) chol2inv(chol(v)))
}

#' @export
bread.glmmTMB <- function(x, ...) {
  cluster <- glmmTMB_components(x)
  X <- model_matrix(x)[, !is.na(coef_CS(x)), drop = FALSE]
  X_list <- matrix_list(X, cluster, "row")
  W_list <- weightMatrix(x, cluster)
  XWX <- Reduce("+", Map(function(X, W) crossprod(X, W %*% X), X_list, W_list))
  chol2inv(chol(XWX)) * v_scale(x)
}
