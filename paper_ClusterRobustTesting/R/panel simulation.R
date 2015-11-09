library(mvtnorm)
library(plyr)
library(tidyr)
library(dplyr)
rm(list = ls())

#------------------------------------------------------
# Set development values for simulation parameters
#------------------------------------------------------

# design matrix parameters

m <- 30 # number of clusters
n <- 20 # number of time-points
k <- 3 # number of outcomes
cluster_balance <- c(A = 1/3, B = 1/3, C = 1/3, AB = 0, AC = 0, BC = 0, ABC = 0) # treatment allocation across clusters
time_balance <- list(A = 1, B = 1, C = 1, AB = c(1/2, 1/2), AC = 0, BC = 0, ABC = c(6 / 10, 3 / 10, 1/10)) # treatment allocation within clusters

# model parameters

rho <- 0.75 # correlation among outcomes (if k > 1)
ar <- 0.6 # auto-correlation across time-points (if k==1)
icc <- 0.3 # cluster-level intra-class correlation
trt_var <- 0.01 # treatment effect variability
outcome_mean <- 0 # average level of the outcome in each treatment condition


#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

within_design <- function(time_balance, n) {
  trt_balance <- sapply(time_balance, function(bal) {
    b <- round(bal * n)
    if (sum(b) > 0) b[1] <- n - sum(b[-1])
    b
  })
  trt_design <- mapply(function(type, balance) rep(type, balance), 
                       type = strsplit(names(time_balance), split = NULL),
                       balance = trt_balance)
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
  if (k == 1) {
    e_mat <- replicate(m, as.vector(arima.sim(list(ar = ar), n = n, 
                                              innov = rnorm(n, sd = sqrt(1 - ar^2)), 
                                              n.start=1, start.innov = rnorm(1))))
    } else {
    Sigma_mat <- rho + diag(1 - rho, nrow = k)
    e_mat <- rmvnorm(m * n, sigma = Sigma_mat)
  }
  return(as.vector(e_mat))
}

simulate_panel <- function(m, n, k, cluster_balance, time_balance, 
                           rho, ar, icc, trt_var, outcome_mean) {

  dat <- design_matrix(m, n, k, cluster_balance, time_balance)
  
  trt_conditions <- unique(dat$trt)
  r <- 1 - trt_var / (2 * icc)
  Tau_mat <- (r + diag(1 - r, nrow = length(trt_conditions))) * icc / (1 - icc)
  
  dat$mu <- unlist(tapply(dat$trt, dat$cluster, function(t) rmvnorm(1, sigma = Tau_mat)[1,t]))
  dat$e <- error_series(m, n, k, rho, ar)
  dat$y <- mu + e
  
  return(dat)
}

select(dat, outcome, cluster, time, trt) %>% spread(time, trt)
select(dat, outcome, cluster, time, y) %>% 
  spread(time, y) %>% 
  group_by(outcome) %>% 
  do(data.frame(r = cor(.[,-1:-2]))) ->
  corr_by_outcome


#------------------------------------------------------
# Model-fitting/estimation/testing functions
#------------------------------------------------------


estimate <- function(dat, design_params) {

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


