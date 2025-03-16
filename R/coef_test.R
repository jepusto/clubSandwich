
#---------------------------------------------
# Satterthwaite approximation
#---------------------------------------------

Satterthwaite_df <- function(P_array) {
  
  V_coef <- 2 * apply(P_array, 3, function(x) sum(x^2))
  E_coef <- apply(P_array, 3, function(x) sum(diag(x)))

  2 * E_coef^2 / V_coef
}

#---------------------------------------------
# Saddlepoint approximation
#---------------------------------------------

saddlepoint_pval <- function(t, Q, eps = 1e-10) {
  if (abs(t) < eps) {
    return(c(s = NA, p_val = 1))
  }
  
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


calc_pval <- function(tstat, df, alternative) {
  switch(
    alternative, 
    `two-sided` = 2 * pt(abs(tstat), df = df, lower.tail = FALSE),
    `greater` = pt(tstat, df = df, lower.tail = FALSE),
    `less` = pt(tstat, df = df, lower.tail = TRUE)
  )
}

#---------------------------------------------
# coeftest for all model coefficients
#---------------------------------------------

#' Test all or selected regression coefficients in a fitted model
#'
#' \code{coef_test} reports one- or two-sided t-tests for each coefficient
#' estimate in a fitted linear regression model, using a sandwich estimator for
#' the standard errors and (optionally) a small sample correction for the
#' p-value. Available small-sample corrections include Satterthwaite
#' approximation or a saddlepoint approximation. Coefficients can be tested
#' against non-zero null values by specifying \code{null_constants}.
#'
#' @param obj Fitted model for which to calculate t-tests.
#' @param vcov Variance covariance matrix estimated using \code{vcovCR} or a
#'   character string specifying which small-sample adjustment should be used to
#'   calculate the variance-covariance.
#' @param test Character vector specifying which small-sample corrections to
#'   calculate. \code{"z"} returns a z test (i.e., using a standard normal
#'   reference distribution). \code{"naive-t"} returns a t test with \code{m -
#'   1} degrees of freedom, where \code{m} is the number of unique clusters.
#'   \code{"naive-tp"} returns a t test with \code{m - p} degrees of freedom,
#'   where \code{p} is the number of regression coefficients in \code{obj}.
#'   \code{"Satterthwaite"} returns a Satterthwaite correction.
#'   \code{"saddlepoint"} returns a saddlepoint correction. Default is
#'   \code{"Satterthwaite"}.
#' @param alternative Character string specifying the alternative hypothesis,
#'   with options "two-sided" (the default), "greater" or "less".
#' @param coefs Character, integer, or logical vector specifying which
#'   coefficients should be tested. The default value \code{"All"} will test all
#'   estimated coefficients.
#' @param null_constants vector of null values for each coefficient to test.
#'   Must have length equal to the number of coefficients specified in
#'   \code{coefs}. Default is \code{0}, in which case the null values are taken
#'   to be zero.
#' @param p_values Logical indicating whether to report p-values. The default
#'   value is \code{TRUE}.
#' @param ... Further arguments passed to \code{\link{vcovCR}}, which are only
#'   needed if \code{vcov} is a character string.
#'
#' @return A data frame containing estimated regression coefficients, standard
#'   errors, specified values of null hypotheses, and test results. For the
#'   Satterthwaite approximation, degrees of freedom and a p-value are reported.
#'   For the saddlepoint approximation, the saddlepoint and a p-value are
#'   reported.
#'
#' @seealso \code{\link{vcovCR}}
#'
#' @examples
#'
#' data("ChickWeight", package = "datasets")
#' lm_fit <- lm(weight ~ Diet  * Time, data = ChickWeight)
#' diet_index <- grepl("Diet.:Time", names(coef(lm_fit)))
#' coef_test(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, coefs = diet_index)
#'
#' V_CR2 <- vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#' coef_test(lm_fit, vcov = V_CR2, coefs = diet_index)
#'
#' # non-inferiority test whether time-by-diet interaction effects are 2 or greater
#' coef_test(lm_fit, vcov = V_CR2, coefs = diet_index, null_constants = 2, alternative = "greater")
#'
#' @export

coef_test <- function(
  obj, 
  vcov, 
  test = "Satterthwaite", 
  alternative = c("two-sided", "greater", "less"), 
  coefs = "All", 
  null_constants = 0, 
  p_values = TRUE, 
  ...
) {
  
  alternative <- match.arg(alternative)
  beta_full <- coef_CS(obj)
  beta_NA <- is.na(beta_full)
  p <- sum(!beta_NA)
  
  which_beta <- get_which_coef(beta_full, coefs)
  
  beta <- beta_full[which_beta & !beta_NA]
  
  if (length(null_constants) == 1L) {
    null_constants <- rep(null_constants, length.out = length(beta))
  }
  if (!is.numeric(null_constants) || length(null_constants) != length(beta)) {
    stop("null_constants must be a numeric vector with length equal to the number of coefficients to be tested.")
  } 
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("z","naive-t","naive-tp","Satterthwaite","saddlepoint")
  if (all(test == "All")) test <- all_tests
  test <- match.arg(test, all_tests, several.ok = TRUE)

  SE <- sqrt(diag(vcov))[which_beta[!beta_NA]]
  
  if (any(c("Satterthwaite","saddlepoint") %in% test)) {
    P_array <- get_P_array(get_GH(obj, vcov))[,,which_beta[!beta_NA],drop=FALSE]
  }
  

  result <- data.frame(Coef = names(beta), beta = as.numeric(beta))
  result$SE <- SE
  result$null_value <- null_constants
  result$tstat <- (beta - null_constants) / SE
  row.names(result) <- result$Coef

  if ("z" %in% test) {
    result$df_z <- Inf
    result$p_z <-  calc_pval(result$tstat, df = Inf, alternative = alternative)
  }
  if ("naive-t" %in% test) {
    J <- nlevels(attr(vcov, "cluster"))
    result$df_t <- J - 1
    result$p_t <-  calc_pval(result$tstat, df = J - 1, alternative = alternative)
  }
  if ("naive-tp" %in% test) {
    J <- nlevels(attr(vcov, "cluster"))
    result$df_tp <- J - p
    result$p_tp <-  calc_pval(result$tstat, df = J - p, alternative = alternative)
  }
  if ("Satterthwaite" %in% test) {
    result$df_Satt <- Satterthwaite_df(P_array = P_array)
    result$p_Satt <- calc_pval(result$tstat, df = result$df_Satt, alternative = alternative)
  }
  if ("saddlepoint" %in% test) {
    saddle <- saddlepoint(t_stats = result$tstat, P_array = P_array)
    result$saddlepoint <- saddle$saddlepoint
    result$p_saddle <- switch(
      alternative,
      `two-sided` = saddle$p_saddle,
      `greater` = ifelse(result$tstat > 0, saddle$p_saddle / 2, 1 - saddle$p_saddle / 2),
      `less` = ifelse(result$tstat > 0, 1 - saddle$p_saddle / 2, saddle$p_saddle / 2)
    )
  }
  
  class(result) <- c("coef_test_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  attr(result, "alternative") <- alternative

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
  
  res$`Null value` <- x$null_value
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
  
  if ("p_tp" %in% names(x)) {
    p_tp <- format.pval(x$p_tp, digits = digits, eps = 10^-digits)
    Sig_tp <- cut(x$p_tp, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res <- cbind(res, "d.f. (naive-tp)" = x$df_tp, "p-val (naive-tp)" = p_tp, "Sig." = Sig_tp)
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

  cat("Alternative hypothesis:", attr(x, "alternative"), "\n")
  print(format(res, digits = 3), row.names = FALSE)
  
}
