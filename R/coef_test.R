
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
# find which coefficients to test
#---------------------------------------------

get_which_coef <- function(beta, coefs) {
  
  p <- length(beta)  
  
  if (identical(coefs,"All")) return(rep(TRUE, p))
  
  switch(class(coefs),
         character = {
           term_names <- names(beta)
           if (length(coefs) == 0) stop("You must specify at least one coefficient to test.")
           if (any(!coefs %in% term_names)) stop("Coefficient names not in model specification.")
           term_names %in% coefs
         },
         logical = {
           if (sum(coefs) == 0) stop("You must specify at least one coefficient to test.")
           if (length(coefs) != p) stop(paste0("Coefficient vector must be of length ",p, "."))
           coefs
         },
         numeric = {
           if (any(!(coefs %in% 1:p))) stop(paste0("Coefficient indices must be less than or equal to ",p,"."))
           if (length(coefs) == 0) stop("You must specify at least one coefficient to test.")
           (1:p) %in% coefs
         },
         integer = {
           if (any(!(coefs %in% 1:p))) stop(paste0("Coefficient indices must be less than or equal to ",p,"."))
           if (length(coefs) == 0) stop("You must specify at least one coefficient to test.")
           (1:p) %in% coefs
         }
         )
}


#---------------------------------------------
# coeftest for all model coefficients
#---------------------------------------------

#' Test all or selected regression coefficients in a fitted model
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
#' @param coefs Character, integer, or logical vector specifying which
#'   coefficients should be tested. The default value \code{"All"} will test all
#'   estimated coefficients.
#' @param p_values Logical indicating whether to report p-values. The default value is \code{TRUE}.
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
#' @examples 
#' data("Produc", package = "plm")
#' lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#' individual_index <- !grepl("state", names(coef(lm_individual)))
#' coef_test(lm_individual, vcov = "CR2", cluster = Produc$state, coefs = individual_index)
#' 
#' V_CR2 <- vcovCR(lm_individual, cluster = Produc$state, type = "CR2")
#' coef_test(lm_individual, vcov = V_CR2, coefs = individual_index)
#' 
#' @export

coef_test <- function(obj, vcov, test = "Satterthwaite", coefs = "All", p_values = TRUE, ...) {
  
  beta_full <- coef_CS(obj)
  beta_NA <- is.na(beta_full)
  
  which_beta <- get_which_coef(beta_full, coefs)
  
  beta <- beta_full[which_beta & !beta_NA]
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("z","naive-t","Satterthwaite","saddlepoint")
  if (all(test == "All")) test <- all_tests
  test <- match.arg(test, all_tests, several.ok = TRUE)

  SE <- sqrt(diag(vcov))[which_beta[!beta_NA]]
  
  if (any(c("Satterthwaite","saddlepoint") %in% test)) {
    P_array <- get_P_array(get_GH(obj, vcov))[,,which_beta[!beta_NA],drop=FALSE]
  }
  

  result <- data.frame(Coef = names(beta), beta = as.numeric(beta))
  result$SE <- SE
  result$tstat <- beta / SE
  row.names(result) <- result$Coef

  if ("z" %in% test) {
    result$df_z <- Inf
    result$p_z <-  2 * pnorm(abs(result$tstat), lower.tail = FALSE)
  }
  if ("naive-t" %in% test) {
    J <- nlevels(attr(vcov, "cluster"))
    result$df_t <- J - 1
    result$p_t <-  2 * pt(abs(result$tstat), df = J - 1, lower.tail = FALSE)
  }
  if ("Satterthwaite" %in% test) {
    Satt <- Satterthwaite(beta = beta, SE = SE, P_array = P_array)
    result$df_Satt <- Satt$df
    result$p_Satt <- Satt$p_Satt
  }
  if ("saddlepoint" %in% test) {
    saddle <- saddlepoint(t_stats = beta / SE, P_array = P_array)
    result$saddlepoint <- saddle$saddlepoint
    result$p_saddle <-saddle$p_saddle
  }
  
  class(result) <- c("coef_test_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")

  if (p_values) {
    result
  } else {
    which_vars <- !grepl("p_", names(result))
    result[which_vars]
  }
  
}

#---------------------------------------------
# print method for coef_test
#---------------------------------------------

#' @export

print.coef_test_clubSandwich <- function(x, digits = 3, ...) {
  
  res <- data.frame(
    `Coef.` = x$Coef, 
    `Estimate` = x$beta, 
    `SE` = x$SE
  )
  
  res$`t-stat` <- x$tstat

  
  if ("p_z" %in% names(x)) {
    p_z <- format.pval(x$p_z, digits = digits, eps = 10^-digits)
    Sig_z <- cut(x$p_z, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res <- cbind(res, "d.f. (z)" = x$df_z,"p-val (z)" = p_z, "Sig." = Sig_z)
  }
  
  if ("p_t" %in% names(x)) {
    p_t <- format.pval(x$p_t, digits = digits, eps = 10^-digits)
    Sig_t <- cut(x$p_t, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res <- cbind(res, "d.f. (naive-t)" = x$df_t, "p-val (naive-t)" = p_t, "Sig." = Sig_t)
  }
  
  if ("p_Satt" %in% names(x)) {
    p_Satt <- format.pval(x$p_Satt, digits = digits, eps = 10^-digits)
    Sig_Satt <- cut(x$p_Satt, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res <- cbind(res, "d.f. (Satt)" = x$df_Satt, "p-val (Satt)" = p_Satt, "Sig." = Sig_Satt)    
  }
  
  if ("p_saddle" %in% names(x)) {
    p_saddle <- format.pval(x$p_saddle, digits = digits, eps = 10^-digits)
    Sig_saddle <- cut(x$p_saddle, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                    labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res <- cbind(res, "s.p." = x$saddlepoint, "p-val (Saddle)" = p_saddle, "Sig." = Sig_saddle)    
  } 

  print(format(res, digits = 3), row.names = FALSE)
  
}
