
#--------------------------------------------------
# confidence intervals for all model coefficients
#---------------------------------------------------

#' Calculate confidence intervals for all or selected regression coefficients in a fitted model
#'
#' \code{conf_int} reports confidence intervals for each coefficient estimate in a fitted
#' linear regression model, using a sandwich estimator for the standard errors
#' and a small sample correction for the critical values. The small-sample correction is
#' based on a Satterthwaite approximation.
#'
#' @param obj Fitted model for which to calculate confidence intervals.
#' @param level Desired coverage level for confidence intervals. 
#' @inheritParams coef_test
#' 
#' @return A data frame containing estimated regression coefficients, standard errors, and confidence intervals. 
#'
#' @seealso \code{\link{vcovCR}}
#' 
#' @examples 
#' data("Produc", package = "plm")
#' lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
#' individual_index <- !grepl("state", names(coef(lm_individual)))
#' conf_int(lm_individual, vcov = "CR2", cluster = Produc$state, coefs = individual_index)
#' 
#' V_CR2 <- vcovCR(lm_individual, cluster = Produc$state, type = "CR2")
#' conf_int(lm_individual, vcov = V_CR2, level = .99, coefs = individual_index)
#'
#' @export

conf_int <- function(obj, vcov, level = .95, test = "Satterthwaite", coefs = "All", ..., p_value = FALSE) {
  
  if (level <= 0 | level >= 1) stop("Confidence level must be between 0 and 1.")
  
  beta_full <- coef_CS(obj)
  beta_NA <- is.na(beta_full)
  
  which_beta <- get_which_coef(beta_full, coefs)
  
  beta <- beta_full[which_beta & !beta_NA]
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("z","naive-t","Satterthwaite")
  test <- match.arg(test, all_tests, several.ok = FALSE)

  SE <- sqrt(diag(vcov))[which_beta[!beta_NA]]

  if (test=="Satterthwaite") {
    P_array <- get_P_array(get_GH(obj, vcov))[,,which_beta[!beta_NA],drop=FALSE]
  }
  
  df <- switch(test, 
               z = Inf,
               `naive-t` = nlevels(attr(vcov, "cluster")) - 1,
               `Satterthwaite` = Satterthwaite(beta = beta, SE = SE, P_array = P_array)$df
  )

  crit <- qt(1 - (1 - level) / 2, df = df)
  
  result <- data.frame(
    beta = beta, 
    SE = SE,
    df = df,
    CI_L = beta - SE * crit,
    CI_U = beta + SE * crit
  )
  
  if (p_value) {
    t_stat <- result$beta / result$SE
    result$p_val <- 2 * pt(abs(t_stat), df = result$df, lower.tail = FALSE)
  }

  class(result) <- c("conf_int_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  attr(result, "level") <- level
  result
}

#---------------------------------------------
# print method for conf_int
#---------------------------------------------

#' @export

print.conf_int_clubSandwich <- function(x, digits = 3, ...) {
  lev <- paste0(100 * attr(x, "level"), "%")
  res <- data.frame("Coef" = rownames(x), x)
  rownames(res) <- NULL
  res_names <- c("Coef", "Estimate", "SE", "d.f.", paste(c("Lower", "Upper"), lev, "CI"))
  if ("p_val" %in% names(x)) {
    res$p_val <- format.pval(x$p_val, digits = digits, eps = 10^-digits)
    res$Sig <- cut(x$p_val, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    res_names <- c(res_names, "p-value", "Sig.")
  }
  names(res) <- res_names
  print(format(res, digits = 3))
}
