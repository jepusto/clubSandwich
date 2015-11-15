library(mvtnorm)
library(plyr)
# library(tidyr)
# library(dplyr)
devtools::load_all() # load clubSandwich

rm(list = ls())

#------------------------------------------------------
# Set development values for simulation parameters
#------------------------------------------------------

# design matrix parameters

# m <- 30 # number of clusters
# n <- 20 # number of time-points
# k <- 3 # number of outcomes
# cluster_balance <- c(A = 1/3, B = 1/3, C = 1/3, AB = 0, AC = 0, BC = 0, ABC = 0) # treatment allocation across clusters
# time_balance <- list(A = 1, B = 1, C = 1, AB = c(1/2, 1/2), AC = 0, BC = 0, ABC = c(6 / 10, 3 / 10, 1/10)) # treatment allocation within clusters
# 
# # model parameters
# 
# rho <- 0.75 # correlation among outcomes (if k > 1)
# ar <- 0.6 # auto-correlation across time-points (if k==1)
# icc <- 0.3 # cluster-level intra-class correlation
# trt_var <- 0.01 # treatment effect variability
# outcome_mean <- rep(0, 3) # average level of the outcome in each treatment condition


#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

within_design <- function(time_balance, n) {
  trt_balance <- lapply(time_balance, function(bal) {
    b <- round(bal * n)
    if (sum(b) > 0) b[1] <- n - sum(b[-1])
    b
  })
  trt_design <- mapply(function(type, balance) rep(type, balance), 
                       type = strsplit(names(time_balance), split = NULL),
                       balance = trt_balance, SIMPLIFY = FALSE)
  names(trt_design) <- names(time_balance)
  trt_design
}

design_matrix <- function(m, n, k, cluster_balance, time_balance) {
  cluster_design <- round(m * cluster_balance)
  cluster_design[1] <- m - sum(cluster_design[-1])
  cluster_type <- rep(names(cluster_design), cluster_design)

  trt_design <- within_design(time_balance, n)
  
  X <- data.frame(outcome = factor(rep(1:k, each = m * n)),
                  cluster = rep(1:m, each = n), 
                  time = 1:n, 
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

simulate_outcome <- function(m, n, k, trt, cluster, icc, rho, ar, trt_var, outcome_mean) {
  
  if (nlevels(trt) != length(outcome_mean)) outcome_mean <- rep_len(outcome_mean, length.out = nlevels(trt))
  
  if (trt_var > 0) {
    r <- 1 - trt_var / (2 * icc)
    Tau_mat <- (r + diag(1 - r, nrow = length(outcome_mean))) * icc / (1 - icc)
    mu <- unlist(tapply(trt, cluster, function(t) rmvnorm(1, mean = outcome_mean, sigma = Tau_mat)[1,t]))
  } else {
    mu <- outcome_mean[trt] + rnorm(m, sd = sqrt(icc / (1 - icc)))[cluster]
  }
  e <- error_series(m, n, k, rho, ar)
  
  as.vector(mu + e)
}

simulate_panel <- function(m, n, k, cluster_balance, time_balance, 
                           icc, rho = 0.8, ar = 0, trt_var = 0, outcome_mean = 0) {

  dat <- design_matrix(m, n, k, cluster_balance, time_balance)
  dat$y <- simulate_outcome(m, n, k, trt = dat$trt, cluster = dat$cluster, 
                            icc, rho, ar, trt_var, outcome_mean)

  return(dat)
}

# dat <- simulate_panel(m = 9, n = 30, k = 3, 
#                       cluster_balance = c(A = 1/3, AB = 1/3, ABC = 1/3), 
#                       time_balance = list(A = 1, AB = c(1/2, 1/2), ABC = c(1/4, 1/4, 1/2)), 
#                       rho = 0.8, ar = 0.0, icc = 0.3, trt_var = 0, outcome_mean = rep(0,3))
# 
# select(dat, outcome, cluster, time, trt) %>% spread(time, trt) %>% head(10)
# 
# select(dat, outcome, cluster, time, y) %>% 
#   spread(time, y) %>% 
#   group_by(outcome) %>% 
#   do(data.frame(r = cor(.[,-1:-2]))) %>%
#   as.data.frame()

#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------

absorb <- function(X_mat, cluster) 
  apply(X_mat, 2, function(x) x - tapply(x, cluster, mean)[cluster])

model_matrix.lm_quick <- function(obj) obj$model_matrix
augmented_model_matrix.lm_quick <- function(obj, cluster, inverse_var) obj$augmented_model_matrix
targetVariance.lm_quick <- function(obj) rep(1, obj$rank + obj$df.residual)
weightMatrix.lm_quick <- function(obj) rep(1, obj$rank + obj$df.residual)

full_fit <- function(dat, constraints) {
  
  cluster_mean <- with(dat, tapply(y, cluster, mean))
  time_mean <- with(dat, tapply(y, time, mean))
  y_absorb <- dat$y - cluster_mean[dat$cluster] - time_mean[dat$time]
  
  S <- model.matrix(~ 0 + factor(time), data = dat)[,-1]
  S_absorb <- absorb(S, dat$cluster)
  R <- model.matrix(~ outcome + trt:outcome, data = dat)[,-1]
  R_absorb_T <- absorb(R, dat$cluster)
  R_absorb <- residuals(lm.fit(S_absorb, R_absorb_T))
  lm_fit <- lm.fit(R_absorb, y_absorb)
  lm_fit$model_matrix <- R_absorb
  lm_fit$augmented_model_matrix <- S_absorb
  class(lm_fit) <- "lm_quick"
  
  C_mats <- lapply(constraints, get_constraint_mat, obj = lm_fit)
  
  V_CR1 <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR1")
  Walds_CR1 <- Wald_test(lm_fit, constraints = constraints, vcov = V_CR1, test = "Naive-F")
  # CR1_adjustments <- lapply(Walds_CR1, function(res) as.data.frame(res[,c("delta","df")]))

  V_CR2 <- vcovCR(lm_fit, cluster = dat$cluster, type = "CR2", inverse_var = TRUE)
  Walds_CR2 <- Wald_test(lm_fit, constraints = constraints, vcov = V_CR2, 
                         test = c("Naive-F","HTA","HTB","HTZ"))
  # CR2_adjustments <- lapply(Walds_CR2, function(res) as.data.frame(res[,c("delta","df")]))

  result <- list(lm_fit = lm_fit, 
                 trt = dat$trt,
                 cluster = dat$cluster,
                 time = dat$time,
                 constraints = C_mats, 
                 CR1_estmats = attr(V_CR1,"estmats"),
                 CR2_estmats = attr(V_CR2,"estmats"),
                 CR1_adjustments = Walds_CR1,
                 CR2_adjustments = Walds_CR2)
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
  
  cluster_mean <- tapply(y, res$cluster, mean)
  time_mean <- tapply(y, res$time, mean)
  y_absorb <- y - cluster_mean[res$cluster] - time_mean[res$time]
  
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
  result <- rbind(Walds_CR1, Walds_CR2)
  rownames(result) <- paste(c(rep("CR1", nrow(res$CR1_adjustments[[1]])),
                              rep("CR2", nrow(res$CR2_adjustments[[1]]))),
                            rownames(result))
  return(result)
}

#---------------------------------
# Test the estimation functions
#---------------------------------

# # generate data
# m = 50
# n = 30
# k = 3
# cluster_balance = c(A = 1/3, AB = 1/3, ABC = 1/3)
# time_balance = list(A = 1, AB = c(1/2, 1/2), ABC = c(1/4, 1/4, 1/2))
# rho = 0.8
# ar = 0.0
# icc = 0.3
# trt_var = 0
# outcome_mean = rep(0, 3)
# constraints <- list(t_B = "outcome1:trtB",
#                     t_C = "outcome1:trtC",
#                     F_1 = c("outcome1:trtB", "outcome1:trtC"),
#                     F_B = c("outcome1:trtB","outcome2:trtB","outcome3:trtB"),
#                     F_C = c("outcome1:trtC","outcome2:trtC","outcome3:trtC"),
#                     F_all = c("outcome1:trtB","outcome2:trtB","outcome3:trtB",
#                               "outcome1:trtC","outcome2:trtC","outcome3:trtC"))
# 
# dat <- simulate_panel(m, n, k, cluster_balance, time_balance, rho, ar, icc, trt_var, outcome_mean)
# 
# # compare full fit and quick fit
# 
# res <- full_fit(dat, constraints)
# quick_res <- quick_fit(res, dat$y)
# cbind(dplyr::bind_rows(res$CR1_adjustments), quick_res[1,])
# identical(dplyr::bind_rows(res$CR1_adjustments)$p_val, as.vector(quick_res[1,]))
# cbind(dplyr::bind_rows(res$CR2_adjustments), as.vector(quick_res[-1,]))
# identical(dplyr::bind_rows(res$CR2_adjustments)$p_val, as.vector(quick_res[-1,]))
# 
# 
# # compare full fit and quick fit on new data
# 
# y <- simulate_outcome(m, n, k, trt = res$trt, cluster = res$cluster, 
#                       icc, rho, ar, trt_var, outcome_mean)
# quick_new <- quick_fit(res, y)
# 
# dat_new <- dat
# dat_new$y <- y
# res_new <- full_fit(dat_new, constraints)
# identical(dplyr::bind_rows(res$CR1_adjustments)[,c("delta","df")], dplyr::bind_rows(res_new$CR1_adjustments)[,c("delta","df")])
# identical(dplyr::bind_rows(res$CR2_adjustments)[,c("delta","df")], dplyr::bind_rows(res_new$CR2_adjustments)[,c("delta","df")])
# 
# cbind(dplyr::bind_rows(res_new$CR1_adjustments), quick_new[1,])
# identical(dplyr::bind_rows(res_new$CR1_adjustments)$p_val, as.vector(quick_new[1,]))
# cbind(rownames(res_new$CR2_adjustments[[1]]), 
#       dplyr::bind_rows(res_new$CR2_adjustments), 
#       quick = as.vector(quick_new[-1,]))
# identical(dplyr::bind_rows(res_new$CR2_adjustments)$p_val, as.vector(quick_new[-1,]))
# 
# # check against lm fit
# 
# lm1 <- lm(y ~ 0 + factor(cluster) + factor(time) + outcome + trt:outcome, data = dat)
# coef(lm1)
# Wald_CR1 <- Wald_test(lm1, constraints = constraints, vcov = "CR1", cluster = dat$cluster, test = "Naive-F")
# dplyr::bind_rows(res$CR1_adjustments)
# dplyr::bind_rows(Wald_CR1)
# all.equal(res$CR1_adjustments, Wald_CR1)
# 
# Wald_CR2 <- Wald_test(lm1, constraints = constraints, vcov = "CR2", cluster = dat$cluster, test = c("Naive-F","HTA","HTB","HTZ"))
# cbind(dplyr::bind_rows(res$CR2_adjustments), dplyr::bind_rows(Wald_CR2))
# all.equal(res$CR2_adjustments, Wald_CR2)

#------------------------------------------------------
# performance calculation
#------------------------------------------------------

calculate_error_rate <- function(a, results) {
  error_rate <- as.data.frame(apply(results < a, 1:2, mean))
  error_rate$alpha <- a
  error_rate$test <- rownames(error_rate)
  error_rate
}

#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

run_sim <- function(iterations, m, n, k, 
                   cluster_balance, time_balance, constraints, 
                   rho, ar, icc, trt_var, outcome_mean, 
                   alpha_levels = c(.005,.01,.05,.10), seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  dat <- simulate_panel(m, n, k, cluster_balance, time_balance, rho, ar, icc, trt_var, outcome_mean)
  initial_fit <- full_fit(dat, constraints)
  
  results <- replicate(iterations, {
    y <- simulate_outcome(m, n, k, trt = dat$trt, cluster = dat$cluster, icc, rho, ar, trt_var, outcome_mean)
    quick_fit(initial_fit, y)
  })

  rejection_rates <- lapply(alpha_levels, calculate_error_rate, results = results)
  do.call(rbind, rejection_rates)
}

# # demonstrate the simulation driver
# 
# cluster_balance <- c(A = 1/3, AB = 1/3, ABC = 1/3)
# time_balance <- list(A = 1, AB = c(1/2, 1/2), ABC = c(1/4, 1/4, 1/2))
# constraints <- list(t_B = "outcome1:trtB",
#                     t_C = "outcome1:trtC",
#                     F_1 = c("outcome1:trtB", "outcome1:trtC"),
#                     F_B = c("outcome1:trtB","outcome2:trtB","outcome3:trtB"),
#                     F_C = c("outcome1:trtC","outcome2:trtC","outcome3:trtC"),
#                     F_all = c("outcome1:trtB","outcome2:trtB","outcome3:trtB",
#                               "outcome1:trtC","outcome2:trtC","outcome3:trtC"))
# 
# system.time(
#   sim_res <- run_sim(iterations = 10, m = 50, n = 18, k = 3, 
#                   cluster_balance, time_balance, constraints,
#                   rho = 0.8, ar = 0, icc = 0.3, trt_var = 0, outcome_mean = rep(0,3))
# )

#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20151114)

# varied design factors

m = c(30,50)
n = c(12,18,30)
icc = c(0, 0.2, 0.4)
trt_var = c(0.01, 0.04)
cluster_balance <- list("all-balanced-within" = c(ABC = 1),
                        "all-unbalanced-within" = c(ABC = 1),
                        "balanced-between" = c(A = 1/3, B = 1/3, C = 1/3),
                        "unbalanced-between" = c(A = .5, B = .3, C = .2),
                        "DD-balanced-within" = c(A = 1/2, ABC = 1/2),
                        "unbalanced-DD-balanced-within" = c(A = 2/3, ABC = 1/3),
                        "DD-unbalanced-within" = c(A = 1/2, ABC = 1/2),
                        "unbalanced-DD-unbalanced-within" = c(A = 2/3, ABC = 1/3))
time_balance

design_factors <- list(m = m, n = n, icc = icc, trt_var = trt_var,
                       cluster_balance = cluster_balance) # combine into a design set
params <- expand.grid(design_factors)

# constant parameters
k = 3 
rho = 0.8
ar = 0
outcome_mean = rep(0,3)
constraints <- list(t_B = "outcome1:trtB",
                    t_C = "outcome1:trtC",
                    F_1 = c("outcome1:trtB", "outcome1:trtC"),
                    F_B = c("outcome1:trtB","outcome2:trtB","outcome3:trtB"),
                    F_C = c("outcome1:trtC","outcome2:trtC","outcome3:trtC"),
                    F_all = c("outcome1:trtB","outcome2:trtB","outcome3:trtB",
                              "outcome1:trtC","outcome2:trtC","outcome3:trtC"))

params$iterations <- 5
params$seed <- round(runif(nrow(params)) * 2^30)

# All look right?
sapply(design_factors, length)
nrow(params)
head(params)



#--------------------------------------------------------
# run simulations in serial
#--------------------------------------------------------
library(plyr)

system.time(results <- mdply(params, .fun = runSim))

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

library(plyr)
library(Pusto)
cluster <- start_parallel(source_obj)

system.time(results <- mdply(params, .fun = run_sim, .parallel = TRUE))

stopCluster(cluster)

save(results, file = "Simulation Results.Rdata")


