#----------------------------------------------
# user-facing vcovCR function
#----------------------------------------------

#' Cluster-robust variance-covariance matrix
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @param obj Fitted model for which to calcualte the variance-covariance matrix
#' @param cluster Expression or vector indicating which observations belong to 
#'   the same cluster. For some classes, the cluster will be detected 
#'   automatically if not specified.
#' @param type Character string specifying which small-sample adjustment should 
#'   be used, with available options \code{"CR0"}, \code{"CR1"}, \code{"CR1S"}, 
#'   \code{"CR2"}, or \code{"CR3"}. See "Details" section of
#'   \code{\link{vcovCR}} for further information.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If a vector, the target matrix is assumed to be 
#'   diagonal. If not specified, \code{vcovCR} will attempt to infer a value.
#' @param inverse_var Optional logical indicating whether the weights used in 
#'   fitting the model are inverse-variance. If not specified, \code{vcovCR} 
#'   will attempt to infer a value.
#' @param form Controls the form of the returned matrix. The default 
#'   \code{"sandwich"} will return the sandwich variance-covariance matrix. 
#'   Alternately, setting \code{form = "meat"} will return only the meat of the 
#'   sandwich and setting \code{form = B}, where \code{B} is a matrix of 
#'   appropriate dimension, will return the sandwich variance-covariance matrix 
#'   calculated using \code{B} as the bread.
#' @param ... Additional arguments available for some classes of objects.
#'   
#' @description This is a generic function, with specific methods defined for 
#'   \code{\link[stats]{lm}}, \code{\link[plm]{plm}}, \code{\link[nlme]{gls}}, 
#'   \code{\link[nlme]{lme}}, \code{\link[robumeta]{robu}}, 
#'   \code{\link[metafor]{rma.uni}}, and \code{\link[metafor]{rma.mv}} objects.
#'   
#' @details Several different small sample corrections are available, which run 
#'   parallel with the "HC" corrections for heteroskedasticity-consistent 
#'   variance estimators, as implemented in \code{\link[sandwich]{vcovHC}}. 
#'   The "CR2" adjustment is recommended (Pustejovsky & Tipton, 2017; 
#'   Imbens & Kolesaar, )
#'   See Pustejovsky and Tipton (2017) and Cameron and Miller (2015) for further 
#'   technical details. Available options include: 
#'   \describe{ 
#'   \item{"CR0"}{is the original form of the sandwich estimator 
#'   (Liang & Zeger, 1986), which does not make any small-sample correction.} 
#'   \item{"CR1"}{multiplies CR0 by \code{m / (m - 1)}, where \code{m} is the
#'   number of clusters.} 
#'   \item{"CR1S"}{multiplies CR0 by \code{(m N) / [(m -
#'   1)(N - p)]}, where \code{m} is the number of clusters, \code{N} is the
#'   total number of observations, and \code{p} is the number of covariates.
#'   Some Stata commands use this correction by default.} 
#'   \item{"CR2"}{is the "bias-reduced linearization" adjustment proposed by 
#'   Bell and McCaffrey (2002) and further developed in Pustejovsky and
#'   Tipton (2017). The adjustment is chosen so that the variance-covariance 
#'   estimator is exactly unbiased under a user-specified working model.}
#'   \item{"CR3"}{approximates the leave-one-cluster-out jackknife variance estimator (Bell & MCCaffery, 2002).}
#'   }
#'   
#' @references 
#' Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28(2), 169-181.
#' 
#' Cameron, A. C., & Miller, D. L. (2015). A Practitioner's Guide to Cluster-Robust Inference. \emph{Journal of Human Resources, 50}(2), 317-372. \doi{10.3368/jhr.50.2.317}
#' 
#' Imbens, G. W., & Kolesar, M. (2016). Robust standard errors in small samples: Some practical advice. \emph{Review of Economics and Statistics, 98}(4), 701-712. \doi{	10.1162/rest_a_00552}
#' 
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using generalized linear models. \emph{Biometrika, 73}(1), 13-22. \doi{10.1093/biomet/73.1.13}
#' 
#' Pustejovsky, J. E. & Tipton, E. (2017). Small sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. \emph{Journal of Business and Economic Statistics}. In Press. \doi{10.1080/07350015.2016.1247004}
#' 
#' @return An object of class \code{c("vcovCR","clubSandwich")}, which consists 
#'   of a matrix of the estimated variance of and covariances between the 
#'   regression coefficient estimates. The matrix has several attributes: 
#'   \describe{ \item{type}{indicates which small-sample adjustment was used} 
#'   \item{cluster}{contains the factor vector that defines independent 
#'   clusters} \item{bread}{contains the bread matrix} \item{v_scale}{constant 
#'   used in scaling the sandwich estimator} \item{est_mats}{contains a list of 
#'   estimating matrices used to calculate the sandwich estimator} 
#'   \item{adjustments}{contains a list of adjustment matrices used to calculate
#'   the sandwich estimator} \item{target}{contains the working 
#'   variance-covariance model used to calculate the adjustment matrices. This 
#'   is needed for calculating small-sample corrections for Wald tests.} }
#'   
#' @seealso \code{\link{vcovCR.lm}}, \code{\link{vcovCR.plm}}, 
#'   \code{\link{vcovCR.gls}}, \code{\link{vcovCR.lme}}, 
#'   \code{\link{vcovCR.robu}}, \code{\link{vcovCR.rma.uni}}, 
#'   \code{\link{vcovCR.rma.mv}}
#' 
#' @examples 
#' 
#' # simulate design with cluster-dependence
#' m <- 8
#' cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
#' n <- length(cluster)
#' X <- matrix(rnorm(3 * n), n, 3)
#' nu <- rnorm(m)[cluster]
#' e <- rnorm(n)
#' y <- X %*% c(.4, .3, -.3) + nu + e
#' dat <- data.frame(y, X, cluster, row = 1:n)
#' 
#' # fit linear model
#' lm_fit <- lm(y ~ X1 + X2 + X3, data = dat)
#' vcov(lm_fit)
#' 
#' # cluster-robust variance estimator with CR2 small-sample correction
#' vcovCR(lm_fit, cluster = dat$cluster, type = "CR2")
#' 
#' # compare small-sample adjustments
#' CR_types <- paste0("CR",c("0","1","1S","2","3"))
#' sapply(CR_types, function(type) 
#'        sqrt(diag(vcovCR(lm_fit, cluster = dat$cluster, type = type))))
#' 
#' @export
#' @import stats

vcovCR <- function(obj, cluster, type, target, inverse_var, form, ...) UseMethod("vcovCR")

#' Cluster-robust variance-covariance matrix
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @rdname vcovCR
#' @export

vcovCR.default <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ...) 
  vcov_CR(obj, cluster, type, target, inverse_var, form)

#---------------------------------------------
# Cluster-robust variance estimator
#---------------------------------------------

adjust_est_mats <- function(type, est_mats, adjustments) {
  switch(type,
         CR0 = est_mats,
         CR1 = lapply(est_mats, function(e) e * adjustments),
         CR1S = lapply(est_mats, function(e) e * adjustments),
         CR2 = Map(function(e, a) e %*% a, e = est_mats, a = adjustments),
         CR3 = Map(function(e, a) e %*% a, e = est_mats, a = adjustments),
         CR4 = Map(function(e, a) a %*% e, e = est_mats, a = adjustments))
}

# uses methods residuals_CS(), model_matrix(), weightMatrix(), 
# targetVariance(), bread(), v_scale()

vcov_CR <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ignore_FE = FALSE) {
  
  cluster <- droplevels(as.factor(cluster))
  
  alias <- is.na(coef_CS(obj))
  X <- model_matrix(obj)
  Xp <- projection_matrix(obj)
  if (any(alias)) {
    X <- X[, !alias, drop = FALSE]
    Xp <- Xp[, !alias, drop = FALSE]
  }  
  
  p <- NCOL(X)
  N <- NROW(X)
  
  if (length(cluster) != N) {
    if (class(na.action(obj)) == "omit") {
      cluster <- droplevels(cluster[-na.action(obj)])
    } else {
      stop("Clustering variable must have length equal to nrow(model_matrix(obj)).")
    }
  } 
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  Xp_list <- matrix_list(Xp, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XpW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = Xp_list, w = W_list)
  
  if (is.null(target)) {
    if (inverse_var) {
      Theta_list <- lapply(W_list, function(w) chol2inv(chol(w)))
    } else {
      Theta_list <- targetVariance(obj, cluster)
    }
  } else {
    if (!is.list(target)) {
      Theta_list <- matrix_list(target, cluster, "both")
    } else {
      Theta_list <- target
    }
  }
  
  if (type %in% c("CR2","CR4")) {
    S <- augmented_model_matrix(obj, cluster, inverse_var, ignore_FE)
    
    if (is.null(S)) {
      rm(S)
      U_list <- Xp_list
      UW_list <- XpW_list
    } else {
      U <- cbind(Xp, S)
      rm(S)
      U_list <- matrix_list(U, cluster, "row")
      UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
    }

    UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
    M_U <- matrix_power(Reduce("+",UWU_list), p = -1)
  }
  
  adjustments <- do.call(type, args = mget(names(formals(type))))
  
  E_list <- adjust_est_mats(type = type, est_mats = XpW_list, adjustments = adjustments)
  
  resid <- residuals_CS(obj)
  res_list <- split(resid, cluster)
  
  components <- do.call(cbind, Map(function(e, r) e %*% r, e = E_list, r = res_list))
  
  v_scale <- v_scale(obj)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) w_scale <- 1L
  
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
  attr(vcov, "est_mats") <- XpW_list
  attr(vcov, "adjustments") <- adjustments
  attr(vcov, "target") <- Theta_list
  attr(vcov, "inverse_var") <- inverse_var
  attr(vcov, "ignore_FE") <- ignore_FE
  class(vcov) <- c("vcovCR","clubSandwich")
  return(vcov)
}

#---------------------------------------------
# as.matrix method for vcovCR
#---------------------------------------------

#' @export

as.matrix.clubSandwich <- function(x, ...) {
  attr(x, "type") <- NULL
  attr(x, "cluster") <- NULL
  attr(x, "bread") <- NULL
  attr(x, "v_scale") <- NULL
  attr(x, "est_mats") <- NULL
  attr(x, "adjustments") <- NULL
  attr(x, "target") <- NULL
  attr(x, "inverse_var") <- NULL
  attr(x, "ignore_FE") <- NULL
  class(x) <- "matrix"
  x
}


#---------------------------------------------
# print method for vcovCR
#---------------------------------------------

#' @export

print.clubSandwich <- function(x, ...) {
  print(as.matrix(x))
}
