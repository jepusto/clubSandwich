library(mvtnorm)
library(plyr)
library(tidyr)
library(dplyr)
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

reps <- 1000
m = 50
n = 30
k = 3
cluster_balance = c(A = 1/3, AB = 1/3, ABC = 1/3)
time_balance = list(A = 1, AB = c(1/2, 1/2), ABC = c(1/4, 1/4, 1/2))
rho = 0.8
ar = 0.0
icc = 0.3
trt_var = 0
outcome_mean = rep(0, 2)
dat <- simulate_panel(m, n, k, cluster_balance, time_balance, rho, ar, icc, trt_var, outcome_mean)
select(dat, outcome, cluster, time, trt) %>% spread(time, trt) %>% head()

select(dat, outcome, cluster, time, y) %>% 
  spread(time, y) %>% 
  group_by(outcome) %>% 
  do(data.frame(r = cor(.[,-1:-2]))) %>%
  as.data.frame()


#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------


full_fit <- function(dat) {
  cluster_means <- with(dat, tapply(y, cluster, mean))
  time_means <- with(dat, tapply(y, time, mean))
  dat$y_absorb <- dat$y - cluster_means[dat$cluster] - time_means[dat$time]
  trt_dummies <- model.matrix(~ trt, data = dat)[,-1]
  dat <- cbind(dat, trt_dummies)
  frml <- paste("y_absorb ~ 0 +", paste(colnames(trt_dummies), "outcome", sep = ":", collapse = " + "))
  lm_fit <- lm(frml, data = dat)
  
  U <- model.matrix(~ 0 + trt + trt:outcome + factor(time), data = dat)
  return(result)
}

re_fit <- function(mod, y) {
  
  return(result)
}

# Test the estimation function

#------------------------------------------------------
# Calculate performance measures
# (For some simulations, it may make more sense
# to do this as part of the simulation driver.)
#------------------------------------------------------

performance <- function(results, model_params) {

  return(performance_measures)
}

# Check performance calculations

#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

runSim <- function(iterations, model_params, design_params, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  results <- replicate(iterations, {
                dat <- dgm(model_params)
                estimate(dat, design_params)
              })

  performance(results, model_params)
}

# demonstrate the simulation driver


#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20150316) # change this seed value!

# now express the simulation parameters as vectors/lists

design_factors <- list(factor1 = , factor2 = , ...) # combine into a design set
params <- expand.grid(design_factors)
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

system.time(results <- mdply(params, .fun = runSim, .parallel = TRUE))

stopCluster(cluster)


save(results, file = "Simulation Results.Rdata")


