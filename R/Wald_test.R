#--------------------------------------------------
# helper functions for constructing constraint matrices
#--------------------------------------------------

#' @name constraint_matrices
#' @title Create constraint matrices
#'
#' @description Helper functions to create common types of constraint matrices,
#'   for use with \code{\link{Wald_test}} to conduct Wald-type tests of linear
#'   contrasts from a fitted regression model.
#'
#' @param constraints Set of constraints to test. Can be logical (using
#'   \code{TRUE} to specify which coefficients to constrain), integer (specify
#'   the index of coefficients to constrain), character (specify the names of
#'   the coefficients to constrain), or a regular expression.
#' @param coefs Vector of coefficient estimates, used to determine the column
#'   dimension of the constraint matrix. Can be omitted if the function is
#'   called inside \code{Wald_test()}.
#' @param reg_ex Logical indicating whether \code{constraints} should be
#'   interpreted as a regular expression. Defaults to \code{FALSE}.
#' @param with_zero Logical indicating whether coefficients should also be
#'   compared to zero. Defaults to \code{FALSE}.
#'
#' @details Constraints can be specified as character vectors, regular
#'   expressions (with \code{reg_ex = TRUE}), integer vectors, or logical
#'   vectors.
#'
#'   \code{constrain_zero()} Creates a matrix that constrains a specified set of
#'   coefficients to all be equal to zero.
#'
#'   \code{constrain_equal()} Creates a matrix that constrains a specified set
#'   of coefficients to all be equal.
#'
#'   \code{constrain_pairwise()} Creates a list of constraint matrices
#'   consisting of all pairwise comparisons between a specified set of
#'   coefficients. If \code{with_zero = TRUE}, then the list will also include a
#'   set of constraint matrices comparing each coefficient to zero.
#'
#' @return A matrix or list of matrices encoding the specified set of
#'   constraints.
#'
#' @seealso \code{\link{Wald_test}}
#'
#' @examples
#' 
#' if (requireNamespace("carData", quietly = TRUE)) withAutoprint({
#' 
#' data(Duncan, package = "carData")
#' Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#'
#' Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
#' # Note that type:income terms are interactions because main effect of income is included
#' # but type:education terms are separate slopes for each unique level of type
#'
#' Duncan_coefs <- coef(Duncan_fit)
#' 
#' # The following are all equivalent
#' constrain_zero(constraints = c("typeprof:income","typewc:income"), 
#'                coefs = Duncan_coefs)
#' constrain_zero(constraints = ":income", coefs = Duncan_coefs, 
#'                reg_ex = TRUE)
#' constrain_zero(constraints = 5:6, coefs = Duncan_coefs)
#' constrain_zero(constraints = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE), 
#'                coefs = Duncan_coefs)
#'
#' # The following are all equivalent
#' constrain_equal(c("typebc:education","typeprof:education","typewc:education"), 
#'                 Duncan_coefs)
#' constrain_equal(":education", Duncan_coefs, reg_ex = TRUE)
#' constrain_equal(7:9, Duncan_coefs)
#' constrain_equal(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE), 
#'                 Duncan_coefs)
#'
#' # Test pairwise equality of the education slopes
#' constrain_pairwise(":education", Duncan_coefs,
#'                    reg_ex = TRUE)
#'
#' # Test pairwise equality of the income slopes, plus compare against zero
#' constrain_pairwise(":income", Duncan_coefs, 
#'                    reg_ex = TRUE, with_zero = TRUE)
#'                    
#' })
#' 
#' @rdname constraint_matrices
#' @export


constrain_zero <- function(constraints, coefs, reg_ex = FALSE) {
  
  if (missing(coefs)) {
    f <- function(coefs) constrain_zero(constraints = constraints, 
                                        coefs = coefs, 
                                        reg_ex = reg_ex)
    return(f)
  }
  
  if (is.list(constraints)) {
    constraint_list <- lapply(constraints, constrain_zero, 
                              coefs = coefs, reg_ex = reg_ex)
    return(constraint_list)
  }
  
  p <- length(coefs)
  
  if (reg_ex) {
    if (!inherits(constraints, "character")) stop("When reg_ex = TRUE, constraints must be a regular expression.")
    constraints <- grepl(constraints, names(coefs))
  }
  
  if ((inherits(constraints, "logical") & sum(as.logical(constraints)) < 1L) | length(constraints) < 1L) stop("You must specify at least one constraint.")

  if (inherits(constraints, "logical")) {
    if (length(constraints) != p) stop(paste0("Constraint logicals must be of length ",p,"."))
    C_mat <- diag(1L, nrow = p)[constraints,,drop=FALSE]
  }
  
  if (inherits(constraints, "numeric") | inherits(constraints, "integer")) {
    if (any(!(constraints %in% 1:p))) stop(paste0("Constraint indices must be less than or equal to ",p,"."))
    C_mat <- diag(1L, nrow = p)[constraints,,drop=FALSE]              
  }
  
  if (inherits(constraints, "character")) {
    term_names <- names(coefs)
    if (any(!constraints %in% term_names)) stop("Constraint names not in model specification.")
    C_mat <- diag(1L, nrow = p)[term_names %in% constraints,,drop=FALSE]
  }

  coef_NA <- is.na(coefs)
  C_mat[,!coef_NA,drop=FALSE]
  
}

#' @rdname constraint_matrices
#' @export

constrain_equal <- function(constraints, coefs, reg_ex = FALSE) {

  if (missing(coefs)) {
    f <- function(coefs) constrain_equal(constraints = constraints, coefs = coefs, reg_ex = reg_ex)
    return(f)
  }
  
  if (is.list(constraints)) {
    constraint_list <- lapply(constraints, constrain_equal, 
                              coefs = coefs, reg_ex = reg_ex)
    return(constraint_list)
  }
  
  if (reg_ex) {
    if (!inherits(constraints, "character")) stop("When reg_ex = TRUE, constraints must be a regular expression.")
    constraints <- grepl(constraints, names(coefs))
  }
  
  if ((inherits(constraints, "logical") & sum(as.logical(constraints)) < 2L) | length(constraints) < 2L) stop("You must specify at least two constraints.")
  
  C_mat <- constrain_zero(constraints = constraints, coefs = coefs)
  
  first_constraint <- which(C_mat[1,] > 0)
  C_mat[,first_constraint] <- -1L
  C_mat[-1,,drop=FALSE]
}

#' @rdname constraint_matrices
#' @export

constrain_pairwise <- function(constraints, coefs, reg_ex = FALSE, with_zero = FALSE) {
  
  if (missing(coefs)) {
    f <- function(coefs) constrain_pairwise(constraints = constraints, 
                                            coefs = coefs, 
                                            reg_ex = reg_ex,
                                            with_zero = with_zero)
    return(f)
  }
  
  if (is.list(constraints)) {
    constraint_list <- lapply(constraints, constrain_pairwise, 
                              coefs = coefs, reg_ex = reg_ex, with_zero = with_zero)
    constraint_list <- unlist(constraint_list, recursive = FALSE)
    return(constraint_list)
  }
  
  p <- length(coefs)
  term_names <- names(coefs)
  
  if (reg_ex) {
    if (!inherits(constraints, "character")) stop("When reg_ex = TRUE, constraints must be a regular expression.")
    constraints <- grepl(constraints, names(coefs))
  }
  
  if ((inherits(constraints, "logical") & sum(as.logical(constraints)) < 2L) | length(constraints) < 2L) stop("You must specify at least two constraints.")
  
  if (inherits(constraints, "logical")) {
    if (length(constraints) != p) stop(paste0("Constraint logicals must be of length ",p,"."))
    constraint_indices <- which(constraints)
  }
  
  if (inherits(constraints, "numeric") | inherits(constraints, "integer")) {
    if (any(!(constraints %in% 1:p))) stop(paste0("Constraint indices must be less than or equal to ",p,"."))
    constraint_indices <- as.integer(constraints)
  }
  
  if (inherits(constraints, "character")) {
    if (!all(constraints %in% term_names)) stop("Constraint names not in model specification.")
    constraint_indices <- which(term_names %in% constraints)
  }
  
  zero_mat <- matrix(0L, nrow = 1, ncol = p)
  
  constraint_pairs <- utils::combn(constraint_indices, 2, simplify = FALSE) 
  
  names(constraint_pairs) <- sapply(constraint_pairs, function(x) paste(term_names[rev(x)], collapse = " - "))
  
  C_mats <- lapply(constraint_pairs, function(x) {
    zero_mat[,x] <- c(-1L, 1L)
    zero_mat
  })
  
  if (with_zero) {
    
    names(constraint_indices) <- term_names[constraint_indices]
    C_to_zero <- lapply(constraint_indices, function(x) {
      zero_mat[,x] <- 1L
      zero_mat
    })
    
    C_mats <- c(C_to_zero, C_mats)
  }
  
  return(C_mats)
  
}


#---------------------------------------------
# Wald-type tests
#---------------------------------------------

#' Test parameter constraints in a fitted linear regression model
#'
#' \code{Wald_test} reports Wald-type tests of linear contrasts from a fitted
#' linear regression model, using a sandwich estimator for the
#' variance-covariance matrix and a small sample correction for the p-value.
#' Several different small-sample corrections are available.
#'
#' @param obj Fitted model for which to calculate Wald tests.
#' @param constraints constraint or list of multiple constraints to test. See
#'   details and examples.
#' @param vcov Variance covariance matrix estimated using \code{vcovCR} or a
#'   character string specifying which small-sample adjustment should be used to
#'   calculate the variance-covariance.
#' @param null_constant vector of null values or list of such vectors for each
#'   set of constraints to test. For a single \code{constraint}, the null values
#'   must have length equal to the number of rows in the constraint. For lists
#'   of null values, each entry must have length equal to the number of rows in
#'   the corresponding entry of \code{constraints}. Default is \code{0}, in
#'   which case the null values are taken to be zero (for every entry, if
#'   \code{constraints} is a list).
#' @param test Character vector specifying which small-sample correction(s) to
#'   calculate. The following corrections are available: \code{"chi-sq"},
#'   \code{"Naive-F"}, \code{"Naive-Fp"}, \code{"HTA"}, \code{"HTB"},
#'   \code{"HTZ"}, \code{"EDF"}, \code{"EDT"}. Default is \code{"HTZ"}.
#' @param tidy Logical value controlling whether to tidy the test results. If
#'   \code{constraints} is a list with multiple constraints, the result will be
#'   coerced into a data frame when \code{tidy = TRUE}.
#' @param ... Further arguments passed to \code{\link{vcovCR}}, which are only
#'   needed if \code{vcov} is a character string.
#'
#' @details Constraints can be specified directly as q X p matrices or
#'   indirectly through \code{\link{constrain_equal}},
#'   \code{\link{constrain_zero}}, or \code{\link{constrain_pairwise}}. By
#'   default, each constraint will be tested against the null hypothesis that it
#'   equal to a zero vector. Non-zero values for null-hypotheses can be
#'   specified using the \code{null_constant} argument.
#'
#' @return A list of test results.
#'
#' @seealso \code{\link{vcovCR}}, \code{\link{constrain_equal}},
#'   \code{\link{constrain_zero}}, \code{\link{constrain_pairwise}}
#'
#' @examples
#'
#'
#' if (requireNamespace("carData", quietly = TRUE)) withAutoprint({
#'
#' data(Duncan, package = "carData")
#' Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#'
#' Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
#' # Note that type:income terms are interactions because main effect of income is included
#' # but type:education terms are separate slopes for each unique level of type
#'
#' # Test equality of intercepts
#' Wald_test(Duncan_fit,
#'           constraints = constrain_equal(1:3),
#'           vcov = "CR2", cluster = Duncan$cluster)
#'
#' # Test equality of type-by-education slopes
#' Wald_test(Duncan_fit,
#'           constraints = constrain_equal(":education", reg_ex = TRUE),
#'           vcov = "CR2", cluster = Duncan$cluster)
#'
#' # Pairwise comparisons of type-by-education slopes
#' Wald_test(Duncan_fit,
#'           constraints = constrain_pairwise(":education", reg_ex = TRUE),
#'           vcov = "CR2", cluster = Duncan$cluster)
#'
#' # Test type-by-income interactions
#' Wald_test(Duncan_fit,
#'           constraints = constrain_zero(":income", reg_ex = TRUE),
#'           vcov = "CR2", cluster = Duncan$cluster)
#'
#' # Pairwise comparisons of type-by-income interactions
#' Wald_test(Duncan_fit,
#'           constraints = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE),
#'           vcov = "CR2", cluster = Duncan$cluster)
#'
#' })
#'
#' @export


Wald_test <- function(obj, constraints, vcov, null_constant = 0, test = "HTZ", tidy = FALSE, ...) {
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!inherits(vcov, "clubSandwich")) stop("Variance-covariance matrix must be a clubSandwich.")
  
  all_tests <- c("chi-sq","Naive-F","Naive-Fp","HTA","HTB","HTZ","EDF","EDT")
  if (all(test == "All")) test <- all_tests
  test <- match.arg(test, all_tests, several.ok = TRUE)

  beta <- na.omit(coef_CS(obj))
  p <- length(beta)
  
  GH <- get_GH(obj, vcov)
  
  # Evaluate constrain_*() functions if used
  if (inherits(constraints, "function")) {
    constraints <- constraints(coef_CS(obj))
  }
  
  if (is.list(constraints)) {
    
    constraints <- lapply(constraints, function(x) {
      if (inherits(x, "function")) x(coef_CS(obj)) else x
    })
    
    # List of constraints
    if (!all(sapply(constraints, inherits, "matrix") & sapply(constraints, ncol) == p)) {
      stop(paste0("Constraints must be a q X ", p," matrix, a list of such matrices, or a call to a constrain_*() function."))
    }
    
    q_constraints <- sapply(constraints, nrow)
    
    null_consts_error_txt <- "Each null_constant must be a numeric vector with length equal to the number of rows in the corresponding constraint matrix."
    
    if (is.list(null_constant)) {
      if (length(null_constant) != length(constraints)) stop("null_constant must be a single vector or a list with the same number of entries as the constraint list.")
      check_numeric <- sapply(null_constant, is.numeric)
      if (!all(check_numeric)) stop(null_consts_error_txt)
    } else if (is.numeric(null_constant)) {
      if (length(null_constant) > 1L) {
        null_constant <- rep(list(null_constant), length(constraints))
      } else {
        null_constant <- lapply(q_constraints, \(x) rep(null_constant, x))
      }
    } else {
      stop(null_consts_error_txt)
    }
    
    q_null <- lengths(null_constant)
    mismatches <- q_constraints != q_null
    if (any(mismatches)) {
      which_mis <- which(mismatches)
      mismatch_msg <- c(
        null_consts_error_txt,
        paste0("Constraint ", which_mis, " has ", q_constraints[mismatches], " rows; null constant ", which_mis, " is length ", q_null[mismatches],".")
      )
      stop(paste(mismatch_msg, collapse = " "))
    }
    
    results <- mapply(
      Wald_testing, 
      C_mat = constraints, null_constant = null_constant, 
      MoreArgs = list(beta = beta, vcov = vcov, test = test, p = p, GH = GH, stop_on_NPD = FALSE),
      SIMPLIFY = FALSE
    )
    
    if (tidy) {
      results <- mapply(
        function(x, nm) cbind(hypothesis = rep(nm, nrow(x)), x, stringsAsFactors = FALSE), 
        x = results, nm = names(results), SIMPLIFY = FALSE
      )
      results <- do.call(rbind, c(results, make.row.names = FALSE))
      class(results) <- c("Wald_test_clubSandwich",class(results))
    }
    
  } else {
    
    
    if (!inherits(constraints, "matrix") | ncol(constraints) != p) {
      stop(paste0("Constraints must be a q X ", p," matrix, a list of such matrices, or a call to a constrain_*() function."))
    }
  
    null_consts_error_txt <- "null_constant must be a numeric vector with length equal to the number of rows in the constraint matrix."
    
    if (is.numeric(null_constant)) {
      if (length(null_constant) == 1L) {
        null_constant <- rep(null_constant, nrow(constraints))
      } else {
        if (length(null_constant) != nrow(constraints)) {
          stop(null_consts_error_txt)
        }
      }
    } else {
      stop(null_consts_error_txt)
    }
    
    results <- Wald_testing(
      C_mat = constraints, null_constant = null_constant, 
      beta = beta, vcov = vcov, test = test, p = p, GH = GH
    ) 
  }
  
  results

}


array_multiply <- function(mat, arr) {
  new_mat <- apply(arr, 3, function(s) mat %*% s)
  array(new_mat, dim = c(nrow(mat), dim(arr)[2], dim(arr)[3]))
}

Wald_testing <- function(C_mat, null_constant, beta, vcov, test, p, GH, stop_on_NPD = TRUE) {
  
  q <- nrow(C_mat)
  dims <- dim(GH$H)
  J <- dims[length(dims)]
  
  if (any(c("HTA","HTB","HTZ","EDF","EDT") %in% test)) {
    GH$G <- lapply(GH$G, function(s) C_mat %*% s)
    if (length(dims)==3) {
      GH$H <- array_multiply(C_mat, GH$H)
    } else {
      H <- array(NA, dim = c(3, q, dims[3:4]))
      for (i in 1:dims[1]) H[i,,,] <- array_multiply(C_mat, GH$H[i,,,])
      GH$H <- H
    }
    P_array <- get_P_array(GH = GH, all_terms = TRUE)
    Omega <- apply(P_array, 1:2, function(x) sum(diag(x)))
    Omega_nsqrt <- matrix_power(Omega, -1/2)
  }
  
  # Wald statistic
  inverse_vcov <- tryCatch(
    chol2inv(chol(C_mat %*% vcov %*% t(C_mat))),
    error = function(e) e
  )
  
  if (inherits(inverse_vcov, "error")) {
    if (stop_on_NPD) {
      stop("Variance-covariance matrix of the contrast is not positive definite. The test cannot be computed.")
    } else {
      result <- data.frame(
        test = test, 
        Fstat = NA_real_, 
        delta = NA_real_, 
        df_num = q, 
        df_denom = NA_real_, 
        p_val = NA_real_
      )
    }
  } else {
    C_beta <- (C_mat %*% beta - matrix(null_constant, ncol = 1L))
    Q <- as.numeric(t(C_beta) %*% inverse_vcov %*% C_beta)
    
    result <- data.frame()
    
    # chi-square
    if ("chi-sq" %in% test) {
      p_val <- pchisq(Q, df = q, lower.tail = FALSE)
      result <- rbind(result, 
                      data.frame(test = "chi-sq", Fstat = Q / q, 
                                 delta = 1, df_num = q, df_denom = Inf, p_val = p_val))
    }
    
    # Naive F
    if ("Naive-F" %in% test) {
      p_val <- pf(Q / q, df1 = q, df2 = J - 1, lower.tail = FALSE)
      result <- rbind(result, 
                      data.frame(test = "Naive-F", Fstat = Q / q, 
                                 delta = 1, df_num = q, df_denom = J - 1, p_val = p_val))
    }
    
    # Naive F with J - p degrees of freedom
    if ("Naive-Fp" %in% test) {
      p_val <- pf(Q / q, df1 = q, df2 = J - p, lower.tail = FALSE)
      result <- rbind(result, 
                      data.frame(test = "Naive-Fp", Fstat = Q / q, 
                                 delta = 1, df_num = q, df_denom = J - p, p_val = p_val))
    }
    
    # Hotelling's T-squared
    if ("HTA" %in% test | "HTB" %in% test) {
      Cov_arr <- covariance_array(P_array, Omega_nsqrt, q = q)
      
      Var_index <- seq(1,q^4, 1 + q^2)
      Var_mat <- matrix(Cov_arr[Var_index], q, q)
      
      
      if ("HTA" %in% test) {
        nu_A <- 2 * sum(Var_mat) / sum(Cov_arr^2)
        result <- rbind(result, data.frame(test = "HTA", Hotelling_Tsq(Q, q, nu = nu_A)))
      } 
      
      if ("HTB" %in% test) {
        lower_mat <- lower.tri(Var_mat, diag = TRUE)
        lower_arr <- array(FALSE, dim = dim(Cov_arr))
        for (s in 1:q) for (t in 1:s) for (u in 1:s) for (v in 1:(ifelse(u==s,t,u))) lower_arr[s,t,u,v] <- TRUE
        
        nu_B <- 2 * sum(Var_mat[lower_mat]) / sum(Cov_arr[lower_arr]^2)
        result <- rbind(result, data.frame(test = "HTB", Hotelling_Tsq(Q, q, nu = nu_B)))
      } 
    } else if ("HTZ" %in% test) {
      Var_mat <- total_variance_mat(P_array, Omega_nsqrt, q = q)
    }
    
    if ("HTZ" %in% test) {
      nu_Z <- q * (q + 1) / sum(Var_mat)
      result <- rbind(result, data.frame(test = "HTZ", Hotelling_Tsq(Q, q, nu = nu_Z)))
    }
    
    # Eigen-decompositions
    
    if ("EDF" %in% test | "EDT" %in% test) {
      spec <- eigen(Omega_nsqrt %*% C_mat %*% vcov %*% t(C_mat) %*% t(Omega_nsqrt))
      df_eig <- 1 / apply(t(spec$vectors) %*% Omega_nsqrt, 1, 
                          function(x) sum(apply(P_array, 3:4, 
                                                function(P) (t(x) %*% P %*% x)^2)))
      
      if ("EDF" %in% test) {
        df4 <- pmax(df_eig, 4.1)
        EQ <- sum(df4 / (df4 - 2))
        VQ <- 2 * sum(df4^2 * (df4 - 1)  / ((df4 - 2)^2 * (df4 - 4))) 
        delta <- ifelse(q * VQ > 2 * EQ^2, (EQ^2 * (q - 2) + 2 * q * VQ) / (EQ * (VQ + EQ^2)), q / EQ)
        df <- ifelse(q * VQ > 2 * EQ^2, 4 + 2 * EQ^2 * (q + 2) / (q * VQ - 2 * EQ^2), Inf)
        Fstat <- delta * Q / q
        p_val <- pf(Fstat, df1 = q, df2 = df, lower.tail = FALSE)
        result <- rbind(result, 
                        data.frame(test = "EDF", Fstat = Fstat, 
                                   delta = delta, df_num = q, df_denom = df, p_val = p_val))
      }
      
      if ("EDT" %in% test) {
        t_j <- t(spec$vectors) %*% Omega_nsqrt %*% C_mat %*% beta / sqrt(spec$values)
        a_j <- df_eig - 1 / 2
        b_j <- 48 * a_j^2
        c_j <- sqrt(a_j * log(1 + t_j^2 / df_eig))
        z_j <- c_j + (c_j^3 + 3 * c_j) / b_j - 
          (4 * c_j^7 + 33 * c_j^5 + 240 * c_j^3 + 855 * c_j) / 
          (10 * b_j^2 + 8 * b_j * c_j^4 + 1000 * b_j)
        Fstat <- mean(z_j^2)
        p_val <- pf(Fstat, df1 = q, df2 = Inf, lower.tail = FALSE)
        result <- rbind(result, 
                        data.frame(test = "EDT", Fstat = Fstat, 
                                   delta = 1, df_num = q, df_denom = Inf, p_val = p_val))
      }
    }
  } 
  
  class(result) <- c("Wald_test_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  result 
}


#--------------------------------------------------
# calculate a covariance array
#--------------------------------------------------

covariance_array <- function(P_array, Omega_nsqrt, q = nrow(Omega_nsqrt)) {
  
  B_jk <- array(apply(P_array, 3:4, function(p) Omega_nsqrt %*% p %*% Omega_nsqrt), 
                dim = dim(P_array))
  
  Cov_arr <- array(NA, dim = rep(q, 4))
  for (s in 1:q) for (t in 1:s) for (u in 1:s) for (v in 1:(ifelse(u==s,t,u))) {
    temp <- sum(B_jk[s,v,,] * B_jk[t,u,,]) + sum(B_jk[s,u,,] * B_jk[t,v,,])
    Cov_arr[s,t,u,v] <- temp
    Cov_arr[s,t,v,u] <- temp
    Cov_arr[t,s,u,v] <- temp
    Cov_arr[t,s,v,u] <- temp
    Cov_arr[u,v,s,t] <- temp
    Cov_arr[u,v,t,s] <- temp
    Cov_arr[v,u,s,t] <- temp
    Cov_arr[v,u,t,s] <- temp
  }
  Cov_arr
}

#---------------------------------------------------------
# calculate total variance of clubSandwich estimator
#---------------------------------------------------------

total_variance_mat <- function(P_array, Omega_nsqrt, q = nrow(Omega_nsqrt)) {
  B_jk <- array(apply(P_array, 3:4, function(p) Omega_nsqrt %*% p %*% Omega_nsqrt), dim = dim(P_array))
  
  var_mat <- matrix(NA, q, q)
  for (s in 1:q) for (t in 1:s) {
    temp <- sum(B_jk[s,t,,] * B_jk[t,s,,]) + sum(B_jk[s,s,,] * B_jk[t,t,,])
    var_mat[s,t] <- temp
    var_mat[t,s] <- temp
  }
  var_mat
}

#--------------------------------------------------
# Hotelling's T-squared approximation
#--------------------------------------------------

Hotelling_Tsq <- function(Q, q, nu) {
  delta <- pmax((nu - q + 1) / nu, 0)
  df <- nu - q + 1
  Fstat <- delta * Q / q
  p_val <- ifelse(df > 0, pf(Fstat, df1 = q, df2 = df, lower.tail = FALSE), as.numeric(NA))
  data.frame(Fstat = Fstat, delta = delta, df_num = q, df_denom = df, p_val = p_val)
}

#---------------------------------------------
# print method for Wald_test
#---------------------------------------------

#' @export

print.Wald_test_clubSandwich <- function(x, digits = 3, ...) {
  res <- x
  res$delta <- NULL
  res$p_val <- format.pval(x$p_val, digits = digits, eps = 10^-digits)
  res$sig <- symnum(x$p_val, corr = FALSE, na = FALSE, 
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                symbols = c("***", "**", "*", ".", " "))
  
  print(format(res, digits = 3), row.names = FALSE)
}
