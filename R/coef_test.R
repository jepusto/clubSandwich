
#---------------------------------------------
# Satterthwaite approximation
#---------------------------------------------

Satterthwaite <- function(beta, SE, S_array, Ex_method = "model") {
  
  V_coef <- 2 * apply(S_array, 1, function(s) sum(crossprod(s)^2))
  
  if (Ex_method == "model") {
    E_coef <- apply(S_array, 1, function(s) sum(s * s))
  } else {
    E_coef <- SE^2
  }
  
  df <- 2 * E_coef^2 / V_coef
  p_val <- 2 * pt(abs(beta / SE), df = df, lower.tail = FALSE)
  data.frame(df = df, p_Satt = p_val)
}

#---------------------------------------------
# Saddlepoint approximation
#---------------------------------------------

saddlepoint_pval <- function(t, Q) {
  eig <- eigen(Q, symmetric = TRUE, only.value=TRUE)$values
  g <- c(1, -t^2 * eig / sum(eig))
  s_eq <- function(s) sum(g / (1 - 2 * g * s))
  s_range <- if (s_eq(0) > 0) c(1 / (2 * min(g)), 0) else c(0, 1 / (2 * max(g)))
  s <- uniroot(s_eq, s_range)$root
  if (s != 0) {
    r <- sign(s) * sqrt(sum(log(1 - 2 * g * s)))
    q <- s * sqrt(2 * sum(g^2 / (1 - 2 * g * s)^2))
    p_val1 <- 1 - pnorm(r) - dnorm(r) * (1 / r - 1 / q)
    p_val2 <- 0.5 - sum(g^3) / (3 * sqrt(pi) * sum(g^2)^(3/2))
    p_val <- min(p_val1, p_val2)
  } else {
    p_val <- 0.5 - sum(g^3) / (3 * sqrt(pi) * sum(g^2)^(3/2))
  }
  c(s = s, p_val = p_val)
}

saddlepoint <- function(t_stats, S_array) {
  saddles <- sapply(1:length(t_stats), function(i) saddlepoint_pval(t = t_stats[i], Q = crossprod(S_array[i,,])))
  data.frame(saddlepoint = saddles["s",], p_saddle = saddles["p_val",])
}

#---------------------------------------------
# coeftest for all model coefficients
#---------------------------------------------

#' Test all regression coefficients in a fitted model
#' 
#' \code{coef_test} reports t-tests for each coefficient estimate in a fitted
#' linear regression model, using a sandwich estimator for the standard errors
#' and a small sample correction for the p-value. The small-sample correction is
#' based on a Satterthwaite approximation or a saddlepoint approximation.
#' 
#' @param obj Fitted model for which to calculate t-tests.
#' @param vcov Variance covariance matrix estimated using \code{vcovCR} or a 
#'   character string specifying which small-sample adjustment should be used to
#'   calculate the variance-covariance.
#' @param test Character vector specifying which small-sample corrections to 
#'   calculate. \code{"Satterthwaite"} returns a Satterthwaite correction. 
#'   \code{"saddlepoint"} returns a saddlepoint correction. Default is 
#'   \code{"Satterthwaite"}.
#' @param Ex_method Character string that controls how the expectation of the 
#'   sandwich estimator is determined for use in the Satterthwaite 
#'   approximation. The default, \code{"model"}, estimates the expectation based
#'   on a working model. The other option, \code{"empirical"} uses the sandwich 
#'   estimate itself.
#' @param ... Further arguments passed to \code{\link{vcovCR}}, which are only 
#'   needed if \code{vcov} is a character string.
#'   
#' @return A data frame containing estimated regression coefficients, standard 
#'   errors, and test results. For the Satterthwaite approximation, degrees of 
#'   freedom and a p-value are reported. For the saddlepoint approximation, the 
#'   saddlepoint and a p-value are reported.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

coef_test <- function(obj, vcov, test = "Satterthwaite", Ex_method = "model", ...) {

  beta <- coef_CR(obj)
  beta_NA <- is.na(beta)
    
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  test <- match.arg(test, c("Satterthwaite","saddlepoint"), several.ok = TRUE)
  Ex_method <- match.arg(Ex_method, c("model","empirical"), several.ok = TRUE)
  
  SE <- sqrt(diag(vcov))
  cluster <- attr(vcov, "cluster")
  E_list <- attr(vcov, "estmats")
  target <- attr(vcov, "target")
  
  S_array <- get_S_array(obj, cluster, target, E_list)
  
  result <- data.frame(beta = beta)
  result$SE[!beta_NA] <- SE
  
  if ("Satterthwaite" %in% test) {
    Satt <- Satterthwaite(beta = beta[!beta_NA], SE = SE, S_array = S_array, Ex_method = Ex_method)
    result$df[!beta_NA] <- Satt$df
    result$p_Satt[!beta_NA] <- Satt$p_Satt
  }
  if ("saddlepoint" %in% test) {
    saddle <- saddlepoint(t_stats = beta[!beta_NA] / SE, S_array = S_array)
    result$saddlepoint[!beta_NA] <- saddle$saddlepoint
    result$p_saddle[!beta_NA] <-saddle$p_saddle
  }
  
  class(result) <- c("coef_test_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  result
}

#---------------------------------------------
# print method for coef_test
#---------------------------------------------

#' @export

print.coef_test_clubSandwich <- function(x, digits = 3, ...) {
  res <- data.frame("Coef" = rownames(x), "Estimate" = x$beta, "Std. Error" = x$SE)
  if ("df" %in% names(x)) {
    p_Satt <- format.pval(x$p_Satt, digits = digits, eps = 10^-digits)

 
    Sig_Satt <- cut(x$p_Satt, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       labels = c("***", "**", "*", ".", " "))
    res <- cbind(res, "d.f." = x$df, "p-val (Satt)" = p_Satt, "Sig." = Sig_Satt)    
  }
  if ("saddlepoint" %in% names(x)) {
    p_saddle <- format.pval(x$p_saddle, digits = digits, eps = 10^-digits)
    Sig_saddle <- cut(x$p_saddle, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    labels = c("***", "**", "*", ".", " "))
    res <- cbind(res, "s.p." = x$saddlepoint, "p-val (Saddle)" = p_saddle, "Sig." = Sig_saddle)    

  } 
  print(format(res, digits = 3))
  res
}

