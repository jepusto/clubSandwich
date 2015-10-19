#--------------------------------------------------
# translate constraint arguments into a matrix
#--------------------------------------------------

get_constraint_mat <- function(obj, constraints) {
  p <- length(coef_CR(obj))
  beta_NA <- is.na(coef_CR(obj))
  
  C_mat <- switch(class(constraints),
        matrix = {
          if (ncol(constraints) != p) stop(paste0("Constraint matrix must have ",p," columns."))
          if (nrow(constraints) == 0) stop("Constraint matrix must have at least one row.")
          constraints       
        },
        logical = {
          if (length(constraints) != p) stop(paste0("Constraint logicals must be of length ",p,"."))
          if (sum(constraints) == 0) stop("You must specify at least one constraint.")
          diag(1L, nrow = p)[constraints,,drop=FALSE]              
        },
        numeric = {
          if (any(!(constraints %in% 1:p))) stop(paste0("Constraint indices must be less than or equal to ",p,"."))
          if (length(constraints) == 0) stop("You must specify at least one constraint.")
          diag(1L, nrow = p)[constraints,,drop=FALSE]              
        },
        integer = {
          if (any(!(constraints %in% 1:p))) stop(paste0("Constraint indices must be less than or equal to ",p,"."))
          if (length(constraints) == 0) stop("You must specify at least one constraint.")
          diag(1L, nrow = p)[constraints,,drop=FALSE]              
        },
        character = {
          term_names <- names(coef_CR(obj))
          if (any(!constraints %in% term_names)) stop("Constraint names not in model specification.")
          if (length(constraints) == 0) stop("You must specify at least one constraint.")
          diag(1L, nrow = p)[term_names %in% constraints,,drop=FALSE]
        })

  C_mat[,!beta_NA,drop=FALSE]

}

#--------------------------------------------------
# calculate a covariance array
#--------------------------------------------------

covariance_array <- function(S_array, Omega_nsqrt, q = nrow(Omega_nsqrt), J = dim(S_array)[3]) {
  B_array <- array(apply(S_array, 3, function(s) Omega_nsqrt %*% s), dim = dim(S_array))
  B_jk <- array(NA, dim = c(J, J, q, q))
  if (q > 1) {
    for (j in 1:J) for (k in 1:j) {
      L <- B_array[,,j] %*% t(B_array[,,k])
      B_jk[j,k,,] <- L
      B_jk[k,j,,] <- t(L)
    }
  } else B_jk[,,1,1] <- crossprod(B_array[1,,])
  
  Cov_arr <- array(NA, dim = rep(q, 4))
  for (s in 1:q) for (t in 1:s) for (u in 1:s) for (v in 1:(ifelse(u==s,t,u))) {
    temp <- sum(B_jk[,,s,v] * B_jk[,,t,u]) + sum(B_jk[,,s,u] * B_jk[,,t,v])
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

total_variance_mat <- function(S_array, Omega_nsqrt, q = nrow(Omega_nsqrt), J = dim(S_array)[3]) {
  B_array <- array(apply(S_array, 3, function(s) Omega_nsqrt %*% s), dim = dim(S_array))
  B_jk <- array(NA, dim = c(J, J, q, q))
  for (j in 1:J) for (k in 1:j) {
    L <- B_array[,,j] %*% t(B_array[,,k])
    B_jk[j,k,,] <- L
    B_jk[k,j,,] <- t(L)
  }
  var_mat <- matrix(NA, q, q)
  for (s in 1:q) for (t in 1:s) {
    temp <- sum(B_jk[,,s,t] * B_jk[,,t,s]) + sum(B_jk[,,s,s] * B_jk[,,t,t])
    var_mat[s,t] <- temp
    var_mat[t,s] <- temp
  }
  var_mat
}

#--------------------------------------------------
# Hotelling's T-squared approximation
#--------------------------------------------------

Hotelling_Tsq <- function(Q, q, nu) {
  delta <- (nu - q + 1) / nu
  df <- nu - q + 1
  Fstat <- delta * Q / q
  p_val <- ifelse(df > 0, pf(Fstat, df1 = q, df2 = df, lower.tail = FALSE), NA)
  c(Fstat = Fstat, delta = delta, df = df, p_val = p_val)
}

#---------------------------------------------
# Wald-type tests
#---------------------------------------------

#' Test all regression coefficients in a fitted model
#' 
#' \code{Wald_test} reports Wald-type tests of linear contrasts from a fitted 
#' linear regression model, using a sandwich estimator for the 
#' variance-covariance matrix and a small sample correction for the p-value. 
#' Several different small-sample corrections are available.
#' 
#' @param obj Fitted model for which to calculate Wald tests.
#' @param constraints List of one or more constraints to test. See details
#'   below.
#' @param vcov Variance covariance matrix estimated using \code{vcovCR} or a 
#'   character string specifying which small-sample adjustment should be used to
#'   calculate the variance-covariance.
#' @param test Character vector specifying which small-sample correction(s) to 
#'   calculate. The following corrections are available: \code{"chi-sq"}, 
#'   \code{"Naive-F"}, \code{"HTA"}, \code{"HTB"}, \code{"HTZ"}, \code{"EDF"}, 
#'   \code{"EDT"}. Default is \code{"HTZ"}.
#' @param ... Further arguments passed to \code{\link{vcovCR}}, which are only 
#'   needed if \code{vcov} is a character string.
#'   
#' @details Constraints can be specified as character vectors, integer vectors,
#' logical vectors, or matrices.
#' 
#' @return A list of test results.
#'   
#' @seealso \code{\link{vcovCR}}
#'   
#' @export

Wald_test <- function(obj, constraints, vcov, test = "HTZ", ...) {
  
  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  if (all(test == "All")) test <- c("chi-sq","Naive-F","HTA","HTB","HTZ","EDF","EDT")
  
  beta <- na.omit(coef_CR(obj))
  
  p <- length(beta)
  
  cluster <- attr(vcov, "cluster")
  E_list <- attr(vcov, "estmats")
  target <- attr(vcov, "target")
  
  S_array <- get_S_array(obj, cluster, target, E_list)
  
  if (is.list(constraints)) {
    C_mats <- lapply(constraints, get_constraint_mat, obj = obj)
    results <- lapply(C_mats, Wald_testing, beta = beta, vcov = vcov, test = test, S_array = S_array)
  } else {
    C_mat <- get_constraint_mat(obj, constraints)
    results <- Wald_testing(C_mat, beta = beta, vcov = vcov, test = test, S_array = S_array) 
  }
  
  results
}

Wald_testing <- function(C_mat, beta, vcov, test, S_array) {
  q <- nrow(C_mat)
  
  if (any(c("chi-sq","Naive-F","HTA","HTB","HTZ","EDF","EDT") %in% test)) {
    N <- dim(S_array)[2]
    J <- dim(S_array)[3]
    S_array <- array(apply(S_array, 3, function(s) C_mat %*% s), dim = c(q, N, J))
    Omega <- apply(array(apply(S_array, 3, tcrossprod), dim = c(q,q,J)), 1:2, sum)
    Omega_nsqrt <- Sym_power(Omega, -1/2)
  }
  
  # Wald statistic
  Q <- as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta)
  
  result <- data.frame(row.names = c("Fstat","delta","df","p_val"))
  
  # chi-square
  if ("chi-sq" %in% test) {
    p_val <- pchisq(Q, df = q, lower.tail = FALSE)
    result <- cbind(result, "chi-sq" = c(Fstat = Q / q, delta = 1, df = Inf, p_val = p_val))
  }
  
  # Naive F
  if ("Naive-F" %in% test) {
    p_val <- pf(Q / q, df1 = q, df2 = J - 1, lower.tail = FALSE)
    result <- cbind(result, "Naive-F" = c(Fstat = Q / q, delta = 1, df = J - 1, p_val = p_val))
  }
  
  # Hotelling's T-squared
  if ("HTA" %in% test | "HTB" %in% test) {
    Cov_arr <- covariance_array(S_array, Omega_nsqrt, q = q, J = J)
    
    Var_index <- seq(1,q^4, 1 + q^2)
    Var_mat <- matrix(Cov_arr[Var_index], q, q)
    
    
    if ("HTA" %in% test) {
      nu_A <- 2 * sum(Var_mat) / sum(Cov_arr^2)
      result <- cbind(result, "HTA" = Hotelling_Tsq(Q, q, nu = nu_A))
    } 
    
    if ("HTB" %in% test) {
      lower_mat <- lower.tri(Var_mat, diag = TRUE)
      lower_arr <- array(FALSE, dim = dim(Cov_arr))
      for (s in 1:q) for (t in 1:s) for (u in 1:s) for (v in 1:(ifelse(u==s,t,u))) lower_arr[s,t,u,v] <- TRUE
      
      nu_B <- 2 * sum(Var_mat[lower_mat]) / sum(Cov_arr[lower_arr]^2)
      result <- cbind(result, "HTB" = Hotelling_Tsq(Q, q, nu = nu_B))
    } 
  } else if ("HTZ" %in% test) {
    Var_mat <- total_variance_mat(S_array, Omega_nsqrt, q = q, J = J)
  }
  
  if ("HTZ" %in% test) {
    nu_Z <- q * (q + 1) / sum(Var_mat)
    result <- cbind(result, "HTZ" = Hotelling_Tsq(Q, q, nu = nu_Z))
  }
  
  # Eigen-decompositions
  
  if ("EDF" %in% test | "EDT" %in% test) {
    spec <- eigen(Omega_nsqrt %*% C_mat %*% vcov %*% t(C_mat) %*% t(Omega_nsqrt))
    D_array <- array(apply(S_array, 3, function(s) t(spec$vector) %*% Omega_nsqrt %*% s), dim = dim(S_array))
    df_eig <- 1 / apply(D_array, 1, function(d) sum(crossprod(d)^2))
    
    if ("EDF" %in% test) {
      df4 <- pmax(df_eig, 4.1)
      EQ <- sum(df4 / (df4 - 2))
      VQ <- 2 * sum(df4^2 * (df4 - 1)  / ((df4 - 2)^2 * (df4 - 4))) 
      delta <- ifelse(q * VQ > 2 * EQ^2, (EQ^2 * (q - 2) + 2 * q * VQ) / (EQ * (VQ + EQ^2)), q / EQ)
      df <- ifelse(q * VQ > 2 * EQ^2, 4 + 2 * EQ^2 * (q + 2) / (q * VQ - 2 * EQ^2), Inf)
      Fstat <- delta * Q / q
      p_val <- pf(Fstat, df1 = q, df2 = df, lower.tail = FALSE)
      result <- cbind(result, "EDF" = c(Fstat = Fstat, delta = delta, df = df, p_val = p_val))
    }
    
    if ("EDT" %in% test) {
      t_j <- t(spec$vector) %*% Omega_nsqrt %*% C_mat %*% beta / sqrt(spec$values)
      a_j <- df_eig - 1 / 2
      b_j <- 48 * a_j^2
      c_j <- sqrt(a_j * log(1 + t_j^2 / df_eig))
      z_j <- c_j + (c_j^3 + 3 * c_j) / b_j - 
        (4 * c_j^7 + 33 * c_j^5 + 240 * c_j^3 + 855 * c_j) / 
        (10 * b_j^2 + 8 * b_j * c_j^4 + 1000 * b_j)
      Fstat <- mean(z_j^2)
      p_val <- pf(Fstat, df1 = q, df2 = Inf, lower.tail = FALSE)
      result <- cbind(result, "EDT" = c(Fstat = Fstat, delta = 1, df = Inf, p_val = p_val))
    }
  }
  
  result <- as.data.frame(t(result))
  class(result) <- c("Wald_test_clubSandwich", class(result))
  attr(result, "type") <- attr(vcov, "type")
  result 
}

#---------------------------------------------
# print method for Wald_test
#---------------------------------------------

#' @export

print.Wald_test_clubSandwich <- function(x, digits = 3, ...) {
  p_val <- format.pval(x$p_val, digits = digits, eps = 10^-digits)
  sig <- symnum(x$p_val, corr = FALSE, na = FALSE, 
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                symbols = c("***", "**", "*", ".", " "))
  res <- data.frame("Test" = rownames(x), "F" = x$Fstat, "d.f." = x$df, "p-val" = p_val)
  #, "Sig." = sig)
  
  print(format(res, digits = 3), row.names = FALSE)
  res
}