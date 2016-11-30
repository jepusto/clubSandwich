
#---------------------------------------------
# Satterthwaite approximation
#---------------------------------------------

Satterthwaite <- function(beta, SE, P_array) {
  
  V_coef <- 2 * apply(P_array, 3, function(x) sum(x^2))
  E_coef <- apply(P_array, 3, function(x) sum(diag(x)))
  
  df <- 2 * E_coef^2 / V_coef
  p_val <- 2 * pt(abs(beta / SE), df = df, lower.tail = FALSE)
  data.frame(df = df, p_Satt = p_val)
}

#---------------------------------------------
# Saddlepoint approximation
#---------------------------------------------

saddlepoint_pval <- function(t, Q) {
  eig <- pmax(0, eigen(Q, symmetric = TRUE, only.values=TRUE)$values)
  g <- c(1, -t^2 * eig / sum(eig))
  s_eq <- function(s) sum(g / (1 - 2 * g * s))
  s_range <- if (t^2 < 1) c(1 / (2 * min(g)), 0) else c(0, 1 / (2 * max(g)))
  s <- uniroot(s_eq, s_range)$root
  if (abs(s) > .01) {
    r <- sign(s) * sqrt(sum(log(1 - 2 * g * s)))
    q <- s * sqrt(2 * sum(g^2 / (1 - 2 * g * s)^2))
    p_val <- 1 - pnorm(r) - dnorm(r) * (1 / r - 1 / q)
  } else {
    p_val <- 0.5 - sum(g^3) / (3 * sqrt(pi) * sum(g^2)^(3/2))
  }
  c(s = s, p_val = p_val)
}

saddlepoint <- function(t_stats, P_array) {
  saddles <- sapply(1:length(t_stats), function(i) saddlepoint_pval(t = t_stats[i], Q = P_array[,,i]))
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
#'   calculate. \code{"z"} returns a z test (i.e., using a standard normal
#'   reference distribution). \code{"naive-t"} returns a t test with \code{m -
#'   1} degrees of freedom. \code{"Satterthwaite"} returns a Satterthwaite
#'   correction. \code{"saddlepoint"} returns a saddlepoint correction. Default
#'   is \code{"Satterthwaite"}.
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

coef_test <- function(obj, vcov, test = "Satterthwaite", ...) {
  
  beta <- coef_CS(obj)
  beta_NA <- is.na(beta)
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("z","naive-t","Satterthwaite","saddlepoint")
  if (all(test == "All")) test <- all_tests
  test <- match.arg(test, all_tests, several.ok = TRUE)
  
  SE <- sqrt(diag(vcov))
  
  if (any(c("Satterthwaite","saddlepoint") %in% test)) {
    P_array <- get_P_array(get_GH(obj, vcov))
  }
  
  result <- data.frame(beta = beta)
  result$SE[!beta_NA] <- SE
  
  if ("z" %in% test) {
    result$p_z[!beta_NA] <-  2 * pnorm(abs(beta[!beta_NA] / SE), lower.tail = FALSE)
  }
  if ("naive-t" %in% test) {
    J <- nlevels(attr(vcov, "cluster"))
    result$p_t[!beta_NA] <-  2 * pt(abs(beta[!beta_NA] / SE), df = J - 1, lower.tail = FALSE)
  }
  if ("Satterthwaite" %in% test) {
    Satt <- Satterthwaite(beta = beta[!beta_NA], SE = SE, P_array = P_array)
    result$df[!beta_NA] <- Satt$df
    result$p_Satt[!beta_NA] <- Satt$p_Satt
  }
  if ("saddlepoint" %in% test) {
    saddle <- saddlepoint(t_stats = beta[!beta_NA] / SE, P_array = P_array)
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
  res <- data.frame("Coef" = rownames(x), "Estimate" = x$beta, "SE" = x$SE)
  if ("p_z" %in% names(x)) {
    p_z <- format.pval(x$p_z, digits = digits, eps = 10^-digits)
    Sig_z <- cut(x$p_z, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels = c("***", "**", "*", ".", " "))
    res <- cbind(res, "p-val (z)" = p_z, "Sig." = Sig_z)
  }
  if ("p_t" %in% names(x)) {
    p_t <- format.pval(x$p_t, digits = digits, eps = 10^-digits)
    Sig_t <- cut(x$p_t, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    labels = c("***", "**", "*", ".", " "))
    res <- cbind(res, "p-val (naive-t)" = p_t, "Sig." = Sig_t)
  }
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
}

