#-----------------------------------------------------
# check that bread can be re-constructed from X and W
#-----------------------------------------------------

check_bread <- function(obj, cluster, y, check_coef = TRUE, tol = 10^-6) {
  cluster <- droplevels(as.factor(cluster))
  B <- sandwich::bread(obj) / v_scale(obj)
  X_list <- matrix_list(model_matrix(obj), cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XWX <- Reduce("+", Map(function(x, w) t(x) %*% w %*% x, x = X_list, w = W_list))
  M <- chol2inv(chol(XWX))
  attr(M, "dimnames") <- attr(B, "dimnames")
  
  eq_bread <- diff(range((B / M)[XWX != 0])) < tol
  
  if (check_coef) {
    coef <- coef_CS(obj)
    y_list <- split(y, cluster)
    XWy <- Reduce("+", Map(function(x, w, y) t(x) %*% w %*% y, x = X_list, w = W_list, y = y_list))
    beta <- as.vector(M %*% XWy)
    names(beta) <- names(coef)
    
    eq_coef <- all.equal(beta, coef, tol = tol)
    if (all(c(eq_coef, eq_bread) == TRUE)) TRUE else list(M = M, B = B, beta = beta, coef = coef)
  } else {
    if (eq_bread) TRUE else list(M = M, B = B)
  }
}

#----------------------------------------------
# check that CR2 and CR4 are target-unbiased
#----------------------------------------------

check_CR <- function(obj, vcov, ..., tol = .Machine$double.eps^0.5) {

  if (is.character(vcov)) vcov <- vcovCR(obj, type = vcov, ...)
  if (!("clubSandwich" %in% class(vcov))) stop("Variance-covariance matrix must be a clubSandwich.")

  # calculate E(V^CRj)  
  cluster <- attr(vcov, "cluster")
  S_array <- get_S_array(obj, vcov)
  E_CRj <- lapply(1:nlevels(cluster), function(j) tcrossprod(S_array[,,j]))
         
  # calculate target
  Theta_list <- attr(vcov, "target")
  X <- model_matrix(obj)
  alias <- is.na(coef_CS(obj))
  if (any(alias)) X <- X[, !alias, drop = FALSE]
  p <- NCOL(X)
  N <- length(cluster)
  J <- nlevels(cluster)
  
  X_list <- matrix_list(X, cluster, "row")
  W_list <- weightMatrix(obj, cluster)
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, w = W_list)
  M <- attr(vcov, "bread") / attr(vcov, "v_scale")
  attr(M, "dimnames") <- NULL

  MXWTWXM <- Map(function(xw, theta) M %*% as.matrix(xw %*% theta %*% t(xw)) %*% M, 
                    xw = XW_list, theta = Theta_list)
  eq <- all.equal(E_CRj, MXWTWXM, tolerance = tol)
  if (all(eq==TRUE)) TRUE else list(E_CRj = E_CRj, target = MXWTWXM)
}


check_sort_order <- function(obj, dat, cluster = NULL,
                             CR_types = paste0("CR",0:3),
                             tol = 10^-6, tol2 = tol, tol3 = tol) {
  
  re_order <- sample(nrow(dat))
  dat_scramble <- dat[re_order,]
  obj_scramble <- update(obj, data = dat_scramble)

  constraints <- utils::combn(length(coef_CS(obj)), 2, simplify = FALSE)
  
  if (is.null(cluster)) {
    CR_fit <- lapply(CR_types, function(x) vcovCR(obj, type = x))
    CR_scramble <- lapply(CR_types, function(x) vcovCR(obj_scramble, type = x))
    test_fit <- lapply(CR_types, function(x) coef_test(obj, vcov = x, test = "All", p_values = FALSE))
    test_scramble <- lapply(CR_types, function(x) coef_test(obj_scramble, vcov = x, test = "All", p_values = FALSE))
    Wald_fit <- Wald_test(obj, constraints = constraints, vcov = "CR2", test = "All")
    Wald_scramble <- Wald_test(obj_scramble, constraints = constraints, vcov = "CR2", test = "All")
    
  } else {
    CR_fit <- lapply(CR_types, function(x) vcovCR(obj, cluster = dat[[cluster]], type = x))
    CR_scramble <- lapply(CR_types, function(x) vcovCR(obj_scramble, cluster = dat_scramble[[cluster]], type = x))
    test_fit <- lapply(CR_types, function(x) coef_test(obj, vcov = x, cluster = dat[[cluster]], test = "All", p_values = FALSE))
    test_scramble <- lapply(CR_types, function(x) coef_test(obj_scramble, vcov = x, cluster = dat_scramble[[cluster]], test = "All", p_values = FALSE))
    Wald_fit <- Wald_test(obj, constraints = constraints, vcov = "CR2",
                          cluster = dat[[cluster]], test = "All")
    Wald_scramble <- Wald_test(obj_scramble, constraints = constraints, vcov = "CR2", 
                               cluster = dat_scramble[[cluster]], test = "All")
    
  }
  
  testthat::expect_equivalent(CR_fit, CR_scramble, tolerance = tol)
  compare_ttests(test_fit, test_scramble, tol = tol2)
  compare_Waldtests(Wald_fit, Wald_scramble, tol = tol3)
}

compare_ttests <- function(a, b, tol = 10^-6) {
  
  if (!inherits(a,"data.frame")) a <- do.call(rbind, a)
  if (!inherits(b,"data.frame")) b <- do.call(rbind, b)
  
  testthat::expect_equal(a$beta, b$beta, tolerance = tol)
  testthat::expect_equal(a$SE, b$SE, tolerance = tol)
  testthat::expect_equal(a$df, b$df, tolerance = tol)
  testthat::expect_equal(a$saddlepoint, b$saddlepoint, tolerance = tol)
}

compare_Waldtests <- function(a, b, tol = 10^-6) {
  
  if (!inherits(a,"data.frame")) a <- do.call(rbind, a)
  if (!inherits(b,"data.frame")) b <- do.call(rbind, b)
  
  testthat::expect_equal(a$Fstat, b$Fstat, tolerance = tol)
  testthat::expect_equal(a$delta, b$delta, tolerance = tol)
  testthat::expect_equal(a$df, b$df, tolerance = tol)
  testthat::expect_equal(a$p_val, b$p_val, tolerance = tol)
  
}