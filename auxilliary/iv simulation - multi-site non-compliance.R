rm(list = ls())

#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

r_site <- function(n, p_trt = 0.5, delta = 0, p_nt = 0.1, 
                   p_at = p_nt, v_uc = 10, r_uy = 0.8, site_id = NULL) {
  
  Trt <- (1:n) %in% (1:round(p_trt * n))
  
  U <- runif(n)
  
  D0 <- rbinom(n, size = 1, prob = qbeta(U, shape1 = v_uc * p_at, shape2 = v_uc * (1 - p_at)))
  D1 <- rbinom(n, size = 1, prob = qbeta(U, shape1 = v_uc * (1 - p_nt) / (1 - p_at), shape2 = v_uc * p_nt / (1 - p_at)))
  D1 <- ifelse(D0 == 1, 1L, D1)
  D <- ((1 - Trt) * D0 + Trt * D1)
  table(D0, D1)
  
  Y0 <- r_uy * qnorm(U) + rnorm(n, mean = 0, sd = sqrt(1 - r_uy^2))
  
  Y <- Y0 + delta * D  
  
  res <- data.frame(Trt = Trt, D = D, Y = Y)
  if (!is.null(site_id)) res$site <- site_id
  
  res
}

sites <- 10
size_mean <- 50
size_sd <- 6
p_trt <- 0.5
delta <- 0
delta_sd <- 0.05
compliance <- 0.5
comp_sd <- 0.0
r_delta_comp <- 0
v_uc <- 5
r_uy <- 0.8


r_multisite_trial <- function(sites, size_mean = 50, size_sd = 0, p_trt = 0.5, 
                              delta = 0, delta_sd = 0, compliance = 0.8, comp_sd = 0.05,
                              r_delta_comp = 0, v_uc = 10, r_uy = 0.8) {
  
  size_p <- 1 - size_sd^2 / (size_mean - 2)
  n_j <- 2 + rbinom(sites, size = round((size_mean - 2) / size_p), prob = size_p)
  
  if (comp_sd > 0) {
    comp_AB <- compliance * (1 - compliance) / comp_sd^2 - 1
    comp_j <- rbeta(sites, shape1 = compliance * comp_AB, shape2 = (1 - compliance) * comp_AB)  
    delta_j <- delta + delta_sd * (r_delta_comp * comp_j / comp_sd + rnorm(sites, mean = 0, sd = sqrt(1 - r_delta_comp^2)))
  } else {
    comp_j <- rep(compliance, sites)
    delta_j <- rnorm(sites, mean = delta, sd = delta_sd)
  }
  
  p_nt <- (1 - comp_j) / 2
  site_dat <- Map(r_site, n = n_j, p_trt = p_trt, delta = delta_j, 
                  p_nt = p_nt, p_at = p_nt, v_uc = v_uc, r_uy = r_uy, 
                  site_id = factor(1:sites))
  
  do.call(rbind, site_dat)
}

dat <- r_multisite_trial(sites = 10, size_mean = 50, size_sd = 6, p_trt = 0.5, 
                         delta = 0, delta_sd = 0.05, compliance = 0.5, comp_sd = 0.0,
                         r_delta_comp = 0, v_uc = 5, r_uy = 0.8)

with(dat, table(Trt, D, site))

#------------------------------------------------------
# Model-fitting function
#------------------------------------------------------

first_stage <- function(dat) {
  first_stage <- lm(D ~ site + site: Trt, data = dat)
  anova(first_stage)[2, "F value"]
}


iv_est <- function(dat) {
  require(AER, quietly = TRUE)
  require(clubSandwich)
  
  lm_fit <- lm(Y ~ site + D, data = dat)
  lm_CR2 <- coef_test(lm_fit, cluster = dat$site, vcov = "CR2", test = c("z", "Satterthwaite"), coefs = "D")
  rownames(lm_CR2) <- NULL
  
  iv_one <- ivreg(Y ~ site + D | site + Trt, data = dat)
  iv_one_CR0 <- coef_test(iv_one, cluster = dat$site, vcov = "CR0", test = c("z", "Satterthwaite"), coefs = "D")
  rownames(iv_one_CR0) <- NULL
  
  iv_one_CR2 <- coef_test(iv_one, cluster = dat$site, vcov = "CR2", test = c("z", "Satterthwaite"), coefs = "D")
  rownames(iv_one_CR2) <- NULL
  
  iv_many <- ivreg(Y ~ site + D | site + site:Trt, data = dat)
  iv_many_CR0 <- coef_test(iv_many, cluster = dat$site, vcov = "CR0", test = c("z", "Satterthwaite"), coefs = "D")
  rownames(iv_many_CR0) <- NULL
  
  iv_many_CR2 <- coef_test(iv_many, cluster = dat$site, vcov = "CR2", test = c("z", "Satterthwaite"), coefs = "D")
  rownames(iv_many_CR2) <- NULL
  
  data.frame(type = c("lm-CR2","iv-one-CR0","iv-one-CR2","iv-many-CR0","iv-many-CR2"), 
             Fstat = Fstat,
             rbind(lm_CR2, iv_one_CR0, iv_one_CR2, iv_many_CR0, iv_many_CR2))
}

iv_est(dat)

#------------------------------------------------------
# Calculate performance measures
#------------------------------------------------------

check_pvals <- function(x, alpha) {
  res <- lapply(alpha, function(a) mean(x < a))
  res <- as.data.frame(res)
  names(res) <- paste("alpha", alpha, sep = "_")
  list(res)
}

check_pvals(x = runif(10000), alpha = c(0.01, 0.05, 0.10))

#-----------------------------------------------------------
# Simulation Driver
#-----------------------------------------------------------

simulate_IV <- function(replicates, sites, size_mean = 40, size_sd = 4,
                        p_trt = 0.5, delta = 0, delta_sd = 0.05,
                        compliance = 0.8, comp_sd = 0.05, r_delta_comp = 0,
                        v_uc = 10, r_uy = 0.8,
                        alpha = c(.01, .05), seed = NULL) {
  require(dplyr) 
  require(tidyr)
  if (!is.null(seed)) set.seed(seed)
  
  reps <- replicate(replicates, {
    dat <- r_multisite_trial(sites = sites, size_mean = size_mean, size_sd = size_sd, p_trt = p_trt, 
                             delta = delta, delta_sd = delta_sd, compliance = compliance, comp_sd = comp_sd,
                             r_delta_comp = r_delta_comp, v_uc = v_uc, r_uy = r_uy)
    iv_est(dat)
  }, simplify = FALSE)
  
  bind_rows(reps) %>%
    group_by(type) %>%
    mutate(V = SE^2) %>%
    summarise(
      E_beta = mean(beta),
      V_beta = var(beta),
      E_var = mean(V),
      V_var = var(V),
      z = check_pvals(p_z, alpha = alpha),
      Satt = check_pvals(p_Satt, alpha = alpha)
    ) %>%
    gather("test","reject", z, Satt) %>%
    unnest(reject)
}

simulate_IV(replicates = 10, sites = 10, size_mean = 50, size_sd = 0, p_trt = 0.5, 
            delta = 0, delta_sd = 0, compliance = 0.8, comp_sd = 0.05,
            r_delta_comp = 0, v_uc = 10, r_uy = 0.8)


#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20171120)

design_factors <- list(
  sites = seq(10, 50, 10), 
  size_mean = 50,
  size_sd = 6, 
  p_trt = 0.5, 
  delta = 0,
  delta_sd = c(0.05, 0.20),
  compliance = seq(.5, .9, .1), 
  comp_sd = c(.00, .05, .10),
  r_delta_comp = 0, 
  v_uc = 5,
  r_uy = 0.8
)

params <- expand.grid(design_factors)
params$replicates <- 5
params$seed <- round(runif(1) * 2^30) + 1:nrow(params)

# All look right?
lengths(design_factors)
nrow(params)
head(params)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

library(Pusto)
clust <- start_parallel(source_obj = source_obj, register = TRUE)

system.time(results <- plyr::mdply(params, .fun = simulate_IV, .parallel = TRUE))

parallel::stopCluster(clust)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "auxilliary/Multisite IV Simulation Results.Rdata")


#--------------------------------------------------------
# Analyze simulation results
#--------------------------------------------------------
library(ggplot2)
rm(list=ls())
load("auxilliary/IV Multisite Simulation Results.Rdata")
