rm(list = ls())

#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

within_design <- function(time_balance, n) {
  trt_balance <- lapply(time_balance, function(bal) {
    b <- round(bal * n)
    if (sum(b) > 0) b[1] <- n - sum(b[-1])
    b
  })
  
  trt_conditions <- str_split(str_extract(names(time_balance), "[[:upper:]]+"), pattern = "")
  trt_design <- mapply(function(type, balance) rep(type, balance), 
                       type = trt_conditions, balance = trt_balance, SIMPLIFY = FALSE)
  names(trt_design) <- names(time_balance)
  trt_design
}

design_matrix <- function(m, n, k, cluster_balance, time_balance, rotating = FALSE) {
  cluster_design <- round(m * cluster_balance)
  cluster_design[1] <- m - sum(cluster_design[-1])
  cluster_type <- rep(names(cluster_design), cluster_design)
  
  trt_design <- within_design(time_balance, n)
  
  if (rotating) {
    offset <- unlist(lapply(table(cluster_type), function(x) 1:x))
    time <- 1:n + rep(offset, each = n)
  } else {
    time <- 1:n
  }
  X <- data.frame(outcome = factor(rep(1:k, each = m * n)),
                  cluster = rep(1:m, each = n), 
                  time = time, 
                  trt = unlist(trt_design[cluster_type]))
  return(X)
}

error_series <- function(m, n, k, rho, ar) {
  if (k == 1 & ar != 0) {
    e_mat <- replicate(m, as.vector(arima.sim(list(ar = ar), n = n, 
                                              innov = rnorm(n, sd = sqrt(1 - ar^2)), 
                                              n.start=1, start.innov = rnorm(1))))
  } else if (k==1 & ar==0) {
    e_mat <- rnorm(m * n)
  } else {
    Sigma_mat <- rho + diag(1 - rho, nrow = k)
    e_mat <- rmvnorm(m * n, sigma = Sigma_mat)
  }
  return(as.vector(e_mat))
}

simulate_outcome <- function(m, n, k, trt, cluster, outcome_id, icc, rho, ar, trt_var, outcome_mean) {
  
  if (nlevels(trt) != length(outcome_mean)) outcome_mean <- rep_len(outcome_mean, length.out = nlevels(trt))
  
  trt <- trt[outcome_id == levels(outcome_id)[1]]
  cluster <- cluster[outcome_id == levels(outcome_id)[1]]
  
  if (trt_var > 0 & icc > 0) {
    r <- 1 - trt_var / (2 * icc)
    Tau_mat <- (r + diag(1 - r, nrow = length(outcome_mean))) * icc / (1 - icc)
    mu <- unlist(tapply(trt, cluster, function(t) rmvnorm(1, mean = outcome_mean, sigma = Tau_mat)[1,t]))
  } else {
    mu <- outcome_mean[trt] + rnorm(m, sd = sqrt(icc / (1 - icc)))[cluster]
  }
  e <- error_series(m, n, k, rho, ar)
  
  as.vector(mu + e)
}

simulate_panel <- function(m, n, k, cluster_balance, time_balance, rotating = FALSE, 
                           icc = 0, rho = 0.8, ar = 0, trt_var = 0, outcome_mean = 0) {
  
  dat <- design_matrix(m, n, k, cluster_balance, time_balance, rotating = rotating)
  dat$y <- simulate_outcome(m, n, k, 
                            trt = dat$trt, cluster = dat$cluster, outcome_id = dat$outcome,
                            icc, rho, ar, trt_var, outcome_mean)
  
  return(dat)
}

# library(corrplot)
# m = 3
# n = 12
# k = 1
# cluster_balance = c(ABC1 = 1/2, ABC2 = 1/2)
# time_balance = list(ABC1 = c(1/3,1/3,1/3), ABC2 = c(1/2,1/4,1/4))
# rho = 0.8
# icc = 0.3
# trt_var = 0.2
# ar = 0
# outcome_mean = 0
# dat <- simulate_panel(m, n, k, cluster_balance, time_balance, 
#                       icc = icc, rho = rho, trt_var = trt_var, 
#                       outcome_mean = outcome_mean)
# trt <- dat$trt
# cluster <- dat$cluster
# outcome_id <- dat$outcome
# 
# more_dat <- replicate(10000, 
#                       simulate_outcome(m, n, k, 
#                                        trt = trt, cluster = cluster, outcome_id = outcome_id,
#                                        icc, rho, ar = 0, trt_var = trt_var, outcome_mean = 0))
# corr_mat <- cor(t(more_dat))
# corrplot(corr_mat, method = "square", cl.lim = c(min(corr_mat), max(corr_mat)))

#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------
absorb <- function(y, x) residuals(lm.fit(x, y))

model_matrix.lm_quick <- function(obj) obj$model_matrix
augmented_model_matrix.lm_quick <- function(obj, cluster, inverse_var) obj$augmented_model_matrix

full_fit <- function(dat, cluster_effects, time_effects, constraints, 
                     tests = c("Naive-F","HTA","HTB","HTZ")) {
  
  R <- model.matrix(~ outcome + trt:outcome, data = dat)
  
  if (cluster_effects & time_effects) {
    Tmat <- model.matrix(~ 0 + factor(cluster), data = dat)
    S <- absorb(y = model.matrix(~ 0 + factor(time), data = dat)[,-1], x = Tmat)
    R_absorb_T <- absorb(y = R[,-1], x = Tmat)
    R_absorb <- absorb(R_absorb_T, S)
    y_absorb_T <- absorb(y = dat$y, x = Tmat)
    y_absorb <- absorb(y_absorb_T, S)
  } else if (time_effects) {
    S <- model.matrix(~ 0 + factor(time), data = dat)
    R_absorb <- absorb(y = R[,-1], x = S)
    y_absorb <- absorb(y = dat$y, x = S)
  } else if (cluster_effects) {
    Tmat <- model.matrix(~ 0 + factor(cluster), data = dat)
    S <- NULL
    R_absorb <- absorb(y = R[,-1], x = Tmat)
    y_absorb <- absorb(y = dat$y, x = Tmat)
  } else {
    S <- NULL
    R_absorb <- R
    y_absorb <- dat$y
  }
  
  lm_fit <- lm.fit(R_absorb, y_absorb)
  lm_fit$model_matrix <- R_absorb
  lm_fit$nobs <- nrow(dat)
  class(lm_fit) <- "lm_quick"
  
  lm_absorb <- lm_fit
  lm_fit$augmented_model_matrix <- S
  
  C_mats <- lapply(constraints, clubSandwich:::get_constraint_mat, obj = lm_fit)
  
  V_CR1 <- vcovCR.default(lm_fit, cluster = dat$cluster, type = "CR1")
  Walds_CR1 <- Wald_test(lm_fit, constraints = constraints, vcov = V_CR1, test = tests)
  
  V_CR2 <- vcovCR.default(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
  Walds_CR2 <- Wald_test(lm_fit, constraints = constraints, vcov = V_CR2, test = tests)
  
  V_CR3 <- vcovCR.default(lm_fit, cluster = dat$cluster, type = "CR3", inverse_var = TRUE)
  Walds_CR3 <- Wald_test(lm_fit, constraints = constraints, vcov = V_CR3, test = tests)
  
  V_CR2A <- vcovCR.default(lm_absorb, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
  Walds_CR2A <- Wald_test(lm_absorb, constraints = constraints, vcov = V_CR2A, test = tests)
  
  result <- list(lm_fit = lm_fit, 
                 trt = dat$trt,
                 cluster = dat$cluster,
                 time = dat$time,
                 cluster_effects = cluster_effects,
                 time_effects = time_effects,
                 constraints = C_mats, 
                 CR1_estmats = attr(V_CR1,"estmats"),
                 CR2_estmats = attr(V_CR2,"estmats"),
                 CR3_estmats = attr(V_CR3,"estmats"),
                 CR2A_estmats = attr(V_CR2A, "estmats"),
                 CR1_adjustments = Walds_CR1,
                 CR2_adjustments = Walds_CR2,
                 CR3_adjustments = Walds_CR3,
                 CR2A_adjustments = Walds_CR2A)
  return(result)
}

vcovCR_quick <- function(resid, cluster, E_list) {
  res_list <- split(resid, cluster)
  components <- do.call(cbind, mapply(function(e, r) e %*% r, e = E_list, r = res_list, SIMPLIFY = FALSE))
  tcrossprod(components)
}

Wald_quick <- function(Q, q, adjustment) {
  F_stat <- adjustment$delta * Q / q
  res <- pf(F_stat, df1 = q, df2 = adjustment$df, lower.tail = FALSE)
  names(res) <- rownames(adjustment)
  res
}

Wald_test_quick <- function(beta, vcov, constraints, adjustments) {
  Q_stats <- lapply(constraints, function(C_mat) 
    as.numeric(t(C_mat %*% beta) %*% chol2inv(chol(C_mat %*% vcov %*% t(C_mat))) %*% C_mat %*% beta))
  q_dim <- lapply(constraints, nrow)
  p_vals <- mapply(Wald_quick, Q = Q_stats, q = q_dim, adjustment = adjustments, SIMPLIFY = FALSE)
  p_mat <- matrix(unlist(p_vals), length(p_vals[[1]]), length(p_vals))
  rownames(p_mat) <- names(p_vals[[1]])
  colnames(p_mat) <- names(p_vals)
  p_mat
}

quick_fit <- function(res, y) {
  
  if (res$cluster_effects & res$time_effects) {
    y_absorb <- absorb(y = y, model.matrix(~ 0 + factor(res$cluster) + factor(res$time)))
  } else if (res$time_effects) {
    y_absorb <- absorb(y = y, model.matrix(~ 0 + factor(res$time)))
  } else if (res$cluster_effects) {
    y_absorb <- absorb(y = y, model.matrix(~ 0 + factor(res$cluster)))
  } else {
    y_absorb <- y
  }
  
  lm_fit <- lm.fit(res$lm_fit$model_matrix, y_absorb)
  
  V_CR1 <- vcovCR_quick(resid = residuals(lm_fit), 
                        cluster = res$cluster, 
                        E_list = res$CR1_estmats)
  
  Walds_CR1 <- Wald_test_quick(beta = coef(lm_fit),
                               vcov = V_CR1,
                               constraints = res$constraints,
                               adjustments = res$CR1_adjustments)
  
  V_CR2 <- vcovCR_quick(resid = residuals(lm_fit), 
                        cluster = res$cluster, 
                        E_list = res$CR2_estmats)
  
  Walds_CR2 <- Wald_test_quick(beta = coef(lm_fit),
                               vcov = V_CR2,
                               constraints = res$constraints,
                               adjustments = res$CR2_adjustments)
  
  V_CR3 <- vcovCR_quick(resid = residuals(lm_fit), 
                        cluster = res$cluster, 
                        E_list = res$CR3_estmats)
  
  Walds_CR3 <- Wald_test_quick(beta = coef(lm_fit),
                               vcov = V_CR3,
                               constraints = res$constraints,
                               adjustments = res$CR3_adjustments)
  
  V_CR2A <- vcovCR_quick(resid = residuals(lm_fit), 
                         cluster = res$cluster, 
                         E_list = res$CR2A_estmats)
  
  Walds_CR2A <- Wald_test_quick(beta = coef(lm_fit),
                                vcov = V_CR2A,
                                constraints = res$constraints,
                                adjustments = res$CR2A_adjustments)
  
  result <- rbind(Walds_CR1, Walds_CR2, Walds_CR3, Walds_CR2A)
  rownames(result) <- paste(c(rep("CR1", nrow(res$CR1_adjustments[[1]])),
                              rep("CR2", nrow(res$CR2_adjustments[[1]])),
                              rep("CR3", nrow(res$CR3_adjustments[[1]])),
                              rep("CR2A", nrow(res$CR3_adjustments[[1]]))),
                            rownames(result))
  return(result)
}

#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

run_sim <- function(iterations, m, n, k, design, constraints, 
                    rho, ar, icc, trt_var, outcome_mean = 0, 
                    alpha_levels = c(.005,.01,.05,.10), seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # get balance parameters
  if (length(design)==1) design <- design[[1]]
  cluster_balance <- design[["cluster_balance"]]
  time_balance <- design[["time_balance"]]
  cluster_effects <- design[["cluster_effects"]]
  time_effects <- design[["time_effects"]]
  rotating <- design$rotating
  if (is.null(rotating)) rotating <- FALSE
  
  # get constraint list
  if (length(constraints) == 1 & is.null(names(constraints))) constraints <- constraints[[1]]
  
  # setup: simulate data and fit model  
  dat <- simulate_panel(m, n, k, cluster_balance, time_balance, rotating,
                        rho, ar, icc, trt_var, outcome_mean)
  initial_fit <- full_fit(dat, cluster_effects, time_effects, constraints)
  
  # replicate based on initial fit 
  results <- replicate(iterations, {
    y <- simulate_outcome(m, n, k, trt = dat$trt, cluster = dat$cluster, 
                          outcome_id = dat$outcome, icc, rho, ar, trt_var, outcome_mean)
    quick_fit(initial_fit, y)
  })
  
  # get degrees of freedom
  CR1_df <- sapply(initial_fit$CR1_adjustments, function(x) x$df)
  CR2_df <- sapply(initial_fit$CR2_adjustments, function(x) x$df)
  CR3_df <- sapply(initial_fit$CR3_adjustments, function(x) x$df)
  CR2A_df <- sapply(initial_fit$CR2A_adjustments, function(x) x$df)
  rownames(CR2_df) <- rownames(initial_fit$CR2_adjustments[[1]])
  df <- rbind(CR1_df, CR2_df, CR3_df, CR2A_df)
  
  # calculate rejection rates
  error_rate <- sapply(alpha_levels, function(a) as.vector(apply(results < a, 1:2, mean)))
  colnames(error_rate) <- paste0("alpha", alpha_levels)
  
  # format output
  data.frame(expand.grid(test = rownames(results), hypothesis = colnames(results)), 
             error_rate, 
             df = as.vector(df))
}

#-------------------------------------
# Experimental Design
#-------------------------------------

source_obj <- ls()

set.seed(20160714)

# balance specifications

designs <- list(
  "RB-balanced" = list(
    cluster_balance = c(ABC = 1),
    time_balance = list(ABC = c(1/2, 1/3, 1/6)),
    cluster_effects = TRUE,
    time_effects = FALSE),
  "RB-unbalanced" = list(
    cluster_balance = c(ABC1 = 1/2, ABC2 = 1/2),
    time_balance = list(ABC1 = c(1/2, 1/3, 1/6), ABC2 = c(1/3, 5/9, 1/9)),
    cluster_effects = TRUE,
    time_effects = FALSE),
  "CR-balanced" = list(
    cluster_balance = c(A = 1/3, B = 1/3, C = 1/3),
    time_balance = list(A = 1, B = 1, C = 1),
    cluster_effects = FALSE,
    time_effects = TRUE),
  "CR-unbalanced" = list(
    cluster_balance = c(A = .5, B = .3, C = .2),
    time_balance = list(A = 1, B = 1, C = 1),
    cluster_effects = FALSE,
    time_effects = TRUE),
  "DD-balanced" = list(
    cluster_balance = c(A = 1/2, ABC = 1/2),
    time_balance = list(A = 1, ABC = c(1/2, 1/3, 1/6)),
    cluster_effects = TRUE,
    time_effects = TRUE),
  "DD-unbalanced" = list(
    cluster_balance = c(A = 2/3, ABC = 1/3),
    time_balance = list(A = 1, ABC = c(1/2, 1/3, 1/6)),
    cluster_effects = TRUE,
    time_effects = TRUE),
  "SD-balanced" = list(
    cluster_balance = c(A = 1/2, ABC = 1/2),
    time_balance = list(A = 1, ABC = c(1/2, 1/3, 1/6)),
    cluster_effects = TRUE,
    time_effects = TRUE,
    rotating = TRUE),
  "SD-unbalanced" = list(
    cluster_balance = c(A = 2/3, ABC = 1/3),
    time_balance = list(A = 1, ABC = c(1/2, 1/3, 1/6)),
    cluster_effects = TRUE,
    time_effects = TRUE, 
    rotating = TRUE)
)

# constraints to test

constraints <- list(t_B = "outcome1:trtB",
                    t_C = "outcome1:trtC",
                    F_1 = c("outcome1:trtB", "outcome1:trtC"),
                    F_B = c("outcome1:trtB","outcome2:trtB","outcome3:trtB"),
                    F_C = c("outcome1:trtC","outcome2:trtC","outcome3:trtC"),
                    F_all = c("outcome1:trtB","outcome2:trtB","outcome3:trtB",
                              "outcome1:trtC","outcome2:trtC","outcome3:trtC"))

# design parameters

design_factors <- list(design = 1:length(designs),
                       iterations = 50,
                       m = c(15, 30, 50), 
                       n = c(6, 18, 30), 
                       icc = c(0.05, 0.15, 0.25), 
                       trt_var = c(0.00, 0.04, 0.09),
                       k = 3,
                       rho = c(0.2, 0.8),
                       ar = 0.0)

params <- expand.grid(design_factors, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
params$design <- designs
params$constraints <- list(constraints)
params$seed <- round(runif(1) * 2^30) + 1:nrow(params)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------
library(Rmpi)
library(snow)
library(foreach)
library(iterators)
library(doSNOW)
library(plyr)

# set up parallel processing
cluster <- getMPIcluster()
registerDoSNOW(cluster)

clusterEvalQ(cluster, library(mvtnorm))
clusterEvalQ(cluster, library(stringr))
clusterEvalQ(cluster, library(clubSandwich))
clusterExport(cluster, source_obj) 

system.time(results <- mdply(params, .fun = run_sim, .parallel = TRUE))

stopCluster(cluster)

results_clean <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
})
head(results_clean, 30)

save(designs, constraints, params, results, file = "Panel simulation results.Rdata")
