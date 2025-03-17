#-------------------------------------
# vcovCR with defaults
#-------------------------------------

#' Cluster-robust variance-covariance matrix for a geeglm object.
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates from an \code{\link[geepack]{geeglm}} object.
#' 
#' @param cluster Expression or vector indicating which observations belong to
#'   the same cluster. Required for \code{geeglm} objects.
#' @param target Optional matrix or vector describing the working
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4}
#'   adjustment matrices. If a vector, the target matrix is assumed to be
#'   diagonal. If not specified, the target is taken to be the estimated variance function.
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
#' if (requireNamespace("geepack", quietly = TRUE)) {
#'
#'   library(geepack)
#'   data(dietox, package = "geepack")
#'   dietox$Cu <- as.factor(dietox$Cu)
#'   mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))
#'   gee1 <- geeglm(mf, data=dietox, id=Pig, family=poisson("identity"), corstr="ar1")
#'   V_CR <- vcovCR(gee1, cluster = dietox$Pig, type = "CR2")
#'   coef_test(gee1, vcov = V_CR, test = "Satterthwaite")
#'   
#' }
#' 
#' @export

vcovCR.geeglm <- function(obj, cluster, type, target = NULL, inverse_var = NULL, form = "sandwich", ...) {
  if (missing(cluster)) {
    cluster <- as.factor(obj$id)
    names(cluster) <- NULL
  } 
  if (is.null(inverse_var)) inverse_var <- is.null(target)
  vcov_CR(obj, cluster = cluster, type = type, 
          target = target, inverse_var = inverse_var, form = form)
}

# coef()
# nobs()

#-----------------------------------------------
# Model matrix
#-----------------------------------------------

#' @export

model_matrix.geeglm <- function(obj) {
  X <- model.matrix(obj)
  eta <- obj$linear.predictors
  dmu_deta <- obj$family$mu.eta
  d <- as.vector(dmu_deta(eta))
  d * X
}

#-------------------------------------
# residuals
#-------------------------------------

#' @export

residuals_CS.geeglm <- function(obj) {
  residuals(obj, type = "response")
}

#-----------------------------------------------
# Get (model-based) working variance matrix 
#-----------------------------------------------

ar1_cor <- function(n, alpha) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  alpha^exponent
}

get_dist <- function(v) {
  mat_dist <- as.matrix(dist(v, diag = TRUE, upper = TRUE))
  mat_dist
}

other_cor <- function(alpha, n = (1 + sqrt(1 + 4 * 2 * length(alpha))) / 2) {
  x <- matrix(1, nrow = n, ncol = n)
  x[lower.tri(x)] <- alpha
  x[upper.tri(x)] <- t(x)[upper.tri(x)]
  x
}

#' @export

targetVariance.geeglm <- function(obj, cluster) {
  
  idvar <- as.factor(obj$id)
  mu <- fitted.values(obj)
  var_fun <- obj$family$variance
  v <- as.numeric(var_fun(mu))
  w <- weights(obj, type = "prior")
  a <- tapply(v / w, idvar, sqrt)
  aa <- lapply(a, tcrossprod)
  
  if (obj$corstr %in% c("independence", "exchangeable", "ar1", "unstructured", "userdefined", "fixed") == F) {
    stop("Working correlation matrix must be a matrix with the following correlation structures: independence, exchangeable, ar1, unstructured, or userdefined")
  } else if (obj$corstr == "ar1") {
    if (is.null(obj$call$waves)) {
      r <- lapply(obj$geese$clusz, ar1_cor, alpha = obj$geese$alpha)
    } else {
      wave <- eval(obj$call$waves, envir = obj$data)
      wave_vec <- split(wave, ceiling(seq_along(wave) / obj$geese$clusz))
      exponent <- lapply(wave_vec, get_dist)
      get_str <- function(alpha, exponent) {
        alpha_str <- alpha^exponent
        alpha_str
      }
      r <- lapply(exponent, get_str, alpha = obj$geese$alpha)
    }
  } else if (obj$corstr == "unstructured") {
    r <- lapply(obj$geese$clusz, other_cor, alpha = as.numeric(obj$geese$alpha))
  } else if (obj$corstr %in% c("userdefined","fixed")) {
    
    formula_env <- attr(obj$formula, ".Environment")
    if (as.character(obj$call$zcor) %in% objects(formula_env)) {
      zcor <- eval(obj$call$zcor, envir = formula_env)
    } else {
      zcor <- eval(obj$call$zcor, envir = parent.frame())
    }
    
    id_cor <- table(idvar)
    id_cor <- rep(names(id_cor), id_cor * (id_cor - 1) / 2)
    
    if (obj$corstr == "userdefined") {
      alpha <- as.numeric(obj$geese$alpha)
      r_vec <- as.numeric(zcor %*% alpha)
      r <- tapply(r_vec, id_cor, other_cor)
    } else if (obj$corstr == "fixed") {
      r <- tapply(zcor, id_cor, other_cor)
    }   
  }
  
  v <- mapply("*", aa, r, SIMPLIFY = FALSE)
  v <- nest_bdiag(v, crosswalk = data.frame(idvar, as.factor(cluster)))
  v
}

#####################
#-------------------------------------
# Get weighting matrix
#-------------------------------------

ar1_cor_inv <- function(n, alpha) {
  if (n == 1) {
    matrix(1)
  } else {
    r_inv <- diag(c(1,rep(1 + alpha^2, n - 2), 1), nrow = n)
    index <- cbind(2:n, 1:(n-1))
    r_inv[index] <- r_inv[index[,2:1]] <- -alpha
    r_inv
  }
}

exch_inv <- function(n, alpha) {
  diag(1 / (1 - alpha), nrow = n) - alpha / ((1 - alpha) * (alpha * (n - 1) + 1))
}

#' @export

weightMatrix.geeglm <- function(obj, cluster) {
  
  idvar <- as.factor(obj$id)
  
  if (obj$corstr %in% c("independence", "exchangeable", "ar1", "unstructured", "userdefined", "fixed") == F) {
    stop("Working correlation matrix must be a matrix with the following correlation structures: independence, exchangeable, ar1, unstructured, or userdefined")
  } else if (obj$corstr %in% c("unstructured","userdefined", "fixed")) {
    
    # Invert the targetVariance for unstructured or user-defined working models
    V_list <- targetVariance.geeglm(obj, idvar)
    W_list <- lapply(V_list, function(v) chol2inv(chol(v)))
  
  } else {
    
    # Otherwise use analytic formulas for inverse of targetVariance
    
    mu <- fitted.values(obj)
    var_fun <- obj$family$variance
    v <- as.numeric(var_fun(mu))
    w <- weights(obj, type = "prior")
    
    if (obj$corstr %in% c("exchangeable","ar1")) {
      a <- tapply(w / v, idvar, sqrt)
      aa <- lapply(a, tcrossprod)
      
      if (obj$corstr == "ar1") {
      
        if (is.null(obj$call$waves)) {
          r_inv <- lapply(obj$geese$clusz, ar1_cor_inv, alpha = obj$geese$alpha)
        } else {
          wave <- eval(obj$call$waves, envir = obj$data)
          wave_vec <- split(wave, ceiling(seq_along(wave) / obj$geese$clusz))
          exponent <- lapply(wave_vec, get_dist)
          get_str <- function(alpha, exponent) {
            alpha_str <- alpha^exponent
            alpha_str
          }
          r <- lapply(exponent, get_str, alpha = obj$geese$alpha)
          r_inv <- lapply(r, function(x) chol2inv(chol(x)))
        }
      } else if (obj$corstr  == "exchangeable") {
        r_inv <- lapply(obj$geese$clusz, exch_inv, alpha = obj$geese$alpha)
      }
      
      W_list <- mapply("*", aa, r_inv, SIMPLIFY = FALSE)
      
    } else {
      W_list <- matrix_list(w / v, idvar, dim = "both")
    }
  }
  
    W_list <- nest_bdiag(W_list, crosswalk = data.frame(idvar, as.factor(cluster)))

  return(W_list)
  
}

#---------------------------------------
# Get bread matrix and scaling constant
#---------------------------------------

#' @export

bread.geeglm <- function(x, ...) {
  cluster <- droplevels(as.factor(x$id))
  X <- model_matrix(x)
  X_list <- matrix_list(X, cluster, "row")
  W_list <- weightMatrix(x, cluster)
  XWX <- Reduce("+", Map(function(x, w) t(x) %*% w %*% x, x = X_list, w = W_list))
  M <- chol2inv(chol(XWX / v_scale(x)))
  rownames(M) <- colnames(M) <- colnames(X)
  M
}

#' @export

v_scale.geeglm <- function(obj) {
  if (substr(obj$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) {
    dispersion <- 1
  } else {
    wres <- as.vector(residuals(obj, "working")) * weights(obj, "working")
    dispersion <- sum(wres^2)/sum(weights(obj, "working"))
  } 
  as.vector(sum(summary(obj)$df[1:2])) * dispersion
}

