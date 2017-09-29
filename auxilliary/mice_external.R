expandvcov <- function(q,u){
	err <- is.na(q)
	return(u)
    ## if (all(!err)) return(u) k <- length(q) v <- names(q) z <- u for (i in 1:ncol(z)){ if (err[i]) { rbind(z[,],NA,z[,]) j
    ## <- j + 1 up <- } j <- j + 1 z[i,] <- u[j,] z[,i] <- u[,j] }
    
    ## z <- matrix(NA, ncol=k, nrow=k, dimnames = list(v,v)) idx <- (is.na()) j <- 0 for (i in 1:k){ if (err[i]) next j <- j
    ## + 1 z[i,] <- u[j,] z[,i] <- u[,j] } return(z)
}

mice.df <- function(m, lambda, dfcom, method) {
	if (is.null(dfcom)) {
		dfcom <- 999999
		warning("Large sample assumed.")
	}
	lambda[lambda < 1e-04] <- 1e-04
	dfold <- (m-1) / lambda^2
	dfobs <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda)
	df <- dfold * dfobs/(dfold + dfobs)
	if (method != "smallsample")
		df <- dfold
	return(df)
}

pool_jp <- function(vcovs, models, n, method = "smallsample") {
  m <- length(models$analyses)
  fa <- models$analyses[[1]] # fa_orig=getfit(models_lm,1)
  analyses <- models$analyses #analyses_orig=getfit(models_lm)
  k <- length(fa$beta) # length(coef(fa_orig))
  names <- rownames(fa) # names(coef(fa_orig))
  qhat <- matrix(NA, nrow = m, ncol = k, dimnames = list(1:m, 
                                                         names))
  u <- array(NA, dim = c(m, k, k), dimnames = list(1:m, names, 
                                                   names))
  for(i in 1:m) {
    fit <- analyses[[i]] #fit2=analyses2[[1]]
    qhat[i, ] <- fit$beta #coef(fit2)
    ui <- vcovs$analyses[[i]] #vcov(fit2)
    ui <- expandvcov(qhat[i, ], ui)
    if (ncol(ui) != ncol(qhat)) 
      stop("Different number of parameters: coef(fit): ", 
           ncol(qhat), ", vcov(fit): ", ncol(ui))
    u[i, , ] <- array(ui, dim = c(1, dim(ui)))
  }
  qbar <- apply(qhat, 2, mean)
  ubar <- apply(u, c(2, 3), mean)
  e <- qhat - matrix(qbar, nrow = m, ncol = k, byrow = TRUE)
  b <- (t(e) %*% e)/(m - 1)
  t <- ubar + (1 + 1/m) * b
  r <- (1 + 1/m) * diag(b/ubar)
  lambda <- (1 + 1/m) * diag(b/t)
  dfcom <- n - (k+1) #df.residual(object)
  df <- mice.df(m, lambda, dfcom, method)
  fmi <- (r + 2/(df + 3))/(r + 1)
  names(r) <- names(df) <- names(fmi) <- names(lambda) <- names
  fit <- list(call = call, call1 = models$call, call2 = models$call1, 
              nmis = models$nmis, m = m, qhat = qhat, u = u, qbar = qbar, 
              ubar = ubar, b = b, t = t, r = r, dfcom = dfcom, df = df, 
              fmi = fmi, lambda = lambda)
  oldClass(fit) <- c("mipo", oldClass(models))
  return(fit)
}




# obj = xyz
# vcov="CR2"
# test = "Satterthwaite"
# cluster = analytic_data$data$Rand.School.ID

# function (obj, vcov, test = "Satterthwaite", ...) 
# {
    # beta <- obj$qbar #beta <- coef_CS(obj)
    # beta_NA <- is.na(beta)
    # if (is.character(vcov)) 
        # vcov <- vcovCR(obj, type = vcov, ...)
    # if (!("clubSandwich" %in% class(vcov))) 
        # stop("Variance-covariance matrix must be a clubSandwich.")
    # all_tests <- c("z", "naive-t", "Satterthwaite", "saddlepoint")
    # if (all(test == "All")) 
        # test <- all_tests
    # test <- match.arg(test, all_tests, several.ok = TRUE)
    # SE <- sqrt(diag(vcov))
    # if (any(c("Satterthwaite", "saddlepoint") %in% test)) {
        # P_array <- get_P_array(get_GH(obj, vcov))
    # }
    # result <- data.frame(beta = beta)
    # result$SE[!beta_NA] <- SE
    # if ("z" %in% test) {
        # result$p_z[!beta_NA] <- 2 * pnorm(abs(beta[!beta_NA]/SE), 
            # lower.tail = FALSE)
    # }
    # if ("naive-t" %in% test) {
        # J <- nlevels(attr(vcov, "cluster"))
        # result$p_t[!beta_NA] <- 2 * pt(abs(beta[!beta_NA]/SE), 
            # df = J - 1, lower.tail = FALSE)
    # }
    # if ("Satterthwaite" %in% test) {
        # Satt <- Satterthwaite(beta = beta[!beta_NA], SE = SE, 
            # P_array = P_array)
        # result$df[!beta_NA] <- Satt$df
        # result$p_Satt[!beta_NA] <- Satt$p_Satt
    # }
    # if ("saddlepoint" %in% test) {
        # saddle <- saddlepoint(t_stats = beta[!beta_NA]/SE, P_array = P_array)
        # result$saddlepoint[!beta_NA] <- saddle$saddlepoint
        # result$p_saddle[!beta_NA] <- saddle$p_saddle
    # }
    # class(result) <- c("coef_test_clubSandwich", class(result))
    # attr(result, "type") <- attr(vcov, "type")
    # result
# }






# vcov_CR <- function(obj, cluster, type, target = NULL, inverse_var = FALSE, form = "sandwich", ignore_FE = FALSE) {
  
  # cluster <- droplevels(as.factor(cluster))
  
  # alias <- is.na(coef_CS(obj))
  # X <- model_matrix(obj)
  # Xp <- projection_matrix(obj)
  # if (any(alias)) {
    # X <- X[, !alias, drop = FALSE]
    # Xp <- Xp[, !alias, drop = FALSE]
  # }  
  
  # p <- NCOL(X)
  # N <- NROW(X)
  
  # if (length(cluster) != N) {
    # if (class(na.action(obj)) == "omit") {
      # cluster <- droplevels(cluster[-na.action(obj)])
    # } else {
      # stop("Clustering variable must have length equal to nrow(model_matrix(obj)).")
    # }
  # } 
  # J <- nlevels(cluster)
  
  # X_list <- matrix_list(X, cluster, "row")
  # Xp_list <- matrix_list(Xp, cluster, "row")
  # W_list <- weightMatrix(obj, cluster)
  # XpW_list <- Map(function(x, w) as.matrix(t(x) %*% w), x = Xp_list, w = W_list)
  
  # if (is.null(target)) {
    # if (inverse_var) {
      # Theta_list <- lapply(W_list, function(w) chol2inv(chol(w)))
    # } else {
      # Theta_list <- targetVariance(obj, cluster)
    # }
  # } else {
    # if (!is.list(target)) {
      # Theta_list <- matrix_list(target, cluster, "both")
    # } else {
      # Theta_list <- target
    # }
  # }
  
  # if (type %in% c("CR2","CR4")) {
    # S <- augmented_model_matrix(obj, cluster, inverse_var, ignore_FE)
    
    # if (is.null(S)) {
      # rm(S)
      # U_list <- Xp_list
      # UW_list <- XpW_list
    # } else {
      # U <- cbind(Xp, S)
      # rm(S)
      # U_list <- matrix_list(U, cluster, "row")
      # UW_list <- Map(function(u, w) as.matrix(t(u) %*% w), u = U_list, w = W_list)
    # }

    # UWU_list <- Map(function(uw, u) uw %*% u, uw = UW_list, u = U_list)
    # M_U <- matrix_power(Reduce("+",UWU_list), p = -1)
  # }
  
  # adjustments <- do.call(type, args = mget(names(formals(type))))
  
  # E_list <- adjust_est_mats(type = type, est_mats = XpW_list, adjustments = adjustments)
  
  # resid <- residuals_CS(obj)
  # res_list <- split(resid, cluster)
  
  # components <- do.call(cbind, Map(function(e, r) e %*% r, e = E_list, r = res_list))
  
  # v_scale <- v_scale(obj)
  # w_scale <- attr(W_list, "w_scale")
  # if (is.null(w_scale)) w_scale <- 1L
  
  # meat <- tcrossprod(components) * w_scale^2 / v_scale
  
  # if (form == "sandwich") {
    # bread <- sandwich::bread(obj)
  # } else if (form == "meat") {
    # bread <- NULL
  # } else if (is.matrix(form)) {
    # bread <- form
    # form <- "sandwich"
  # } 
  
  # vcov <- switch(form, 
                 # sandwich = bread %*% meat %*% bread / v_scale,
                 # meat = meat)
  # rownames(vcov) <- colnames(vcov) <- colnames(X)
  # attr(vcov, "type") <- type
  # attr(vcov, "cluster") <- cluster
  # attr(vcov, "bread") <- bread
  # attr(vcov, "v_scale") <- v_scale
  # attr(vcov, "est_mats") <- XpW_list
  # attr(vcov, "adjustments") <- adjustments
  # attr(vcov, "target") <- Theta_list
  # attr(vcov, "inverse_var") <- inverse_var
  # attr(vcov, "ignore_FE") <- ignore_FE
  # class(vcov) <- c("vcovCR","clubSandwich")
  # return(vcov)
# }


