
#--------------------------------------------------
# confidence intervals for all model coefficients
#---------------------------------------------------

#' Calculate confidence intervals for all or selected regression coefficients in
#' a fitted model
#'
#' \code{conf_int} reports confidence intervals for each coefficient estimate in
#' a fitted linear regression model, using a sandwich estimator for the standard
#' errors and a small sample correction for the critical values. The
#' small-sample correction is based on a Satterthwaite approximation.
#'
#' @param obj Fitted model for which to calculate confidence intervals.
#' @param level Desired coverage level for confidence intervals.
#' @param test Character vector specifying which small-sample corrections to
#'   calculate. \code{"z"} returns a z test (i.e., using a standard normal
#'   reference distribution). \code{"naive-t"} returns a t test with \code{m -
#'   1} degrees of freedom, where \code{m} is the number of unique clusters.
#'   \code{"naive-tp"} returns a t test with \code{m - p} degrees of freedom,
#'   where \code{p} is the number of regression coefficients in \code{obj}.
#'   \code{"Satterthwaite"} returns a Satterthwaite correction. Unlike in
#'   \code{coef_test()}, \code{"saddlepoint"} is not currently supported in
#'   \code{conf_int()} because saddlepoint confidence intervals do not have a
#'   closed-form solution.
#' @param p_values Logical indicating whether to report p-values. The default
#'   value is \code{FALSE}.

#' @inheritParams coef_test
#'
#' @return A data frame containing estimated regression coefficients, standard
#'   errors, confidence intervals, and (optionally) p-values.
#'
#' @seealso \code{\link{vcovCR}}
#'
#' @examples
#' 
#' data("ChickWeight", package = "datasets")
#' lm_fit <- lm(weight ~ Diet  * Time, data = ChickWeight)
#' diet_index <- grepl("Diet.:Time", names(coef(lm_fit)))
#' conf_int(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, coefs = diet_index)
#'
#' V_CR2 <- vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#' conf_int(lm_fit, vcov = V_CR2, level = .99, coefs = diet_index)
#'
#' @export

conf_int <- function(obj, vcov, level = .95, test = "Satterthwaite", coefs = "All", ..., p_values = FALSE) {
  
  if (level <= 0 | level >= 1) stop("Confidence level must be between 0 and 1.")
  
  beta_full <- coef_CS(obj)
  beta_NA <- is.na(beta_full)
  p <- sum(!beta_NA)
    
  which_beta <- get_which_coef(beta_full, coefs)
  
  beta <- beta_full[which_beta & !beta_NA]
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("z","naive-t","naive-tp","Satterthwaite")
  if (test == "saddlepoint") stop("test = 'saddlepoint' is not currently supported  because saddlepoint confidence intervals do not have a closed-form solution.")
  test <- match.arg(test, all_tests, several.ok = FALSE)

  SE <- sqrt(diag(vcov))[which_beta[!beta_NA]]

  if (test=="Satterthwaite") {
    P_array <- get_P_array(get_GH(obj, vcov))[,,which_beta[!beta_NA],drop=FALSE]
  }
  
  df <- switch(test, 
               z = Inf,
               `naive-t` = nlevels(attr(vcov, "cluster")) - 1,
               `naive-tp` = nlevels(attr(vcov, "cluster")) - p,
               `Satterthwaite` = Satterthwaite(beta = beta, SE = SE, P_array = P_array)$df
  )

  crit <- qt(1 - (1 - level) / 2, df = df)
  
  result <- data.frame(
    Coef = names(beta),
    beta = beta, 
    SE = SE,
    df = df,
    CI_L = beta - SE * crit,
    CI_U = beta + SE * crit
   )
  row.names(result) <- result$Coef

  if (p_values) {
    t_stat <- result$beta / result$SE
    result$p_val <- 2 * pt(abs(t_stat), df = result$df, lower.tail = FALSE)
  }

  class(result) <- c("conf_int_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  attr(result, "level") <- level
  result
}

#--------------------------------------------------
# confidence intervals for linear contrasts
#---------------------------------------------------

#' Calculate confidence intervals and p-values for linear contrasts of
#' regression coefficients in a fitted model
#'
#' \code{linear_contrast} reports confidence intervals and (optionally) p-values
#' for linear contrasts of regression coefficients from a fitted model, using a
#' sandwich estimator for the standard errors and (optionally) a small sample
#' correction for the critical values. The default small-sample correction is
#' based on a Satterthwaite approximation.
#'
#' @param obj Fitted model for which to calculate confidence intervals.
#' @param contrasts A contrast matrix, or a list of multiple contrast matrices
#'   to test. See details and examples.
#' @param level Desired coverage level for confidence intervals.
#' @param test Character vector specifying which small-sample corrections to
#'   calculate. \code{"z"} returns a z test (i.e., using a standard normal
#'   reference distribution). \code{"naive-t"} returns a t test with \code{m -
#'   1} degrees of freedom, where \code{m} is the number of unique clusters.
#'   \code{"naive-tp"} returns a t test with \code{m - p} degrees of freedom,
#'   where \code{p} is the number of regression coefficients in \code{obj}.
#'   \code{"Satterthwaite"} returns a Satterthwaite correction. Unlike in
#'   \code{coef_test()}, \code{"saddlepoint"} is not currently supported in
#'   \code{conf_int()} because saddlepoint confidence intervals do not have a
#'   closed-form solution.
#' @param p_values Logical indicating whether to report p-values. The default
#'   value is \code{FALSE}.
#' @inheritParams coef_test
#' 
#' @details Constraints can be specified directly as q X p matrices or
#'   indirectly through \code{\link{constrain_pairwise}},
#'   \code{\link{constrain_equal}}, or \code{\link{constrain_zero}}.
#'
#' @return A data frame containing estimated contrasts, standard
#'   errors, confidence intervals, and (optionally) p-values.
#'
#' @seealso \code{\link{vcovCR}}
#'
#' @examples
#' 
#' data("ChickWeight", package = "datasets")
#' lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)
#' 
#' # Pairwise comparisons of diet-by-time slopes
#' linear_contrast(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, 
#'                 contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE))
#'
#' 
#' if (requireNamespace("carData", quietly = TRUE)) withAutoprint({
#' 
#'   data(Duncan, package = "carData")
#'   Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#'
#'   Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
#'   # Note that type:income terms are interactions because main effect of income is included
#'   # but type:education terms are separate slopes for each unique level of type
#'
#'   # Pairwise comparisons of type-by-education slopes
#'   linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
#'                   contrasts = constrain_pairwise(":education", reg_ex = TRUE),
#'                   test = "Satterthwaite")
#'  
#'   # Pairwise comparisons of type-by-income interactions
#'   linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
#'                   contrasts = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE),
#'                   test = "Satterthwaite")
#'                   
#' })
#'
#' @export

linear_contrast <- function(obj, vcov, contrasts, level = .95, test = "Satterthwaite", ..., p_values = FALSE) {
  
  if (level <= 0 | level >= 1) stop("Confidence level must be between 0 and 1.")
  
  beta_full <- coef_CS(obj)
  beta_NA <- is.na(beta_full)
  p <- sum(!beta_NA)
  beta <- beta_full[!beta_NA]

  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")

  # Evaluate constrain_*() functions if used
  if (inherits(contrasts, "function")) {
    contrasts <- contrasts(beta_full)
  }
  
  if (is.list(contrasts)) {
    
    contrast_list <- lapply(contrasts, function(x) {
      if (inherits(x, "function")) x(beta_full) else x
    })
    
    if (any(sapply(contrast_list, is.list))) {
      contrast_list <- lapply(contrast_list, function(x) if (is.list(x)) x else list(x))
      contrast_list <- do.call(c, args = contrast_list)
    }

    # List of contrasts
    if (!all(sapply(contrast_list, inherits, "matrix") & sapply(contrast_list, ncol) == p)) {
      stop(paste0("Contrasts must be a q X ", p," matrix, a list of such matrices, or a call to a constrain_*() function."))
    }
    
    # Add row names to contrasts
    contrast_list <- mapply(function(x,r) {
      if (is.null(rownames(x))) {
        if (nrow(x) > 1) {
          rownames(x) <- paste(r, 1:nrow(x), sep = ".")
        } else {
          rownames(x) <- r
        }
      } else {
        rownames(x) <- paste(r, rownames(x), sep = ".")
      }
      return(x) 
      }, 
      x = contrast_list, r = names(contrast_list), 
      SIMPLIFY = FALSE)
    
    # Combine into one matrix
    contrasts <- do.call(rbind, args = contrast_list)
  } else if (is.matrix(contrasts)) {
    if (is.null(rownames(contrasts))) {
      rownames(contrasts) <- paste("Contrast", 1:NROW(contrasts))
    }
  } else if (is.vector(contrasts)) {
    contrasts <- matrix(contrasts, nrow = 1L)
    rownames(contrasts) <- "Contrast"
  }
    
  all_tests <- c("z","naive-t","naive-tp","Satterthwaite")
  if (test == "saddlepoint") stop("test = 'saddlepoint' is not currently supported  because saddlepoint confidence intervals do not have a closed-form solution.")
  test <- match.arg(test, all_tests, several.ok = FALSE)
  
  # Calculate contrasts
  est <- contrasts %*% beta
  
  # Standard error of the contrasts
  vcov_contrasts <- contrasts %*% vcov %*% t(contrasts)
  SE <- sqrt(diag(vcov_contrasts))
  
  # Satterthwaite degrees of freedom
  if (test=="Satterthwaite") {
    q <- nrow(contrasts)
    GH <- get_GH(obj, vcov)
    GH$G <- lapply(GH$G, function(s) contrasts %*% s)
    dims <- dim(GH$H)
    if (length(dims)==3) {
      GH$H <- array_multiply(contrasts, GH$H)
    } else {
      H <- array(NA, dim = c(3, q, dims[3:4]))
      for (i in 1:dims[1]) H[i,,,] <- array_multiply(contrasts, GH$H[i,,,,drop=FALSE])
      GH$H <- H
    }
    P_array <- get_P_array(GH = GH)
  }
  
  df <- switch(test, 
               z = Inf,
               `naive-t` = nlevels(attr(vcov, "cluster")) - 1,
               `naive-tp` = nlevels(attr(vcov, "cluster")) - p,
               `Satterthwaite` = Satterthwaite(beta = est, SE = SE, P_array = P_array)$df
  )
  
  crit <- qt(1 - (1 - level) / 2, df = df)
  
  result <- data.frame(
    Coef = rownames(contrasts),
    Est = est, 
    SE = SE,
    df = df,
    CI_L = est - SE * crit,
    CI_U = est + SE * crit
  )
  row.names(result) <- result$Coef
  
  if (p_values) {
    t_stat <- result$Est / result$SE
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
  res_names <- c("Coef.", "Estimate", "SE", "d.f.", paste(c("Lower", "Upper"), lev, "CI"))
  if ("p_val" %in% names(x)) {
    x$Sig <- cut(x$p_val, breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels = c("***", "**", "*", ".", " "), include.lowest = TRUE)
    x$p_val <- format.pval(x$p_val, digits = digits, eps = 10^-digits)
    res_names <- c(res_names, "p-value", "Sig.")
  }
  names(x) <- res_names
  print(format(x, digits = 3), row.names = FALSE)
}
