rm(list = ls())

#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

# clusters <- 40
# p_trt <- 0.5
# delta <- 0
# size_mean <- 5
# p_nt <- 0.1
# p_at <- 0.1
# v_uc <- 10
# r_uy <- 0.8
# icc <- 0.2

r_cluster_trial <- function(clusters, p_trt = 0.5, 
                            delta = 0, size_mean = 10, 
                            p_nt = 0.1, p_at = p_nt, v_uc = 10, 
                            r_uy = 0.8, icc = 0.2) {
  
  n_j <- 1 + rpois(clusters, lambda = size_mean)
  cluster_id <- rep(1:clusters, n_j)
  Trt <- (1:clusters) %in% (1:round(p_trt * clusters))
  
  U <- runif(clusters)
  D0 <- rbinom(clusters, size = 1, prob = qbeta(U, shape1 = v_uc * p_at, shape2 = v_uc * (1 - p_at)))
  D1 <- rbinom(clusters, size = 1, prob = qbeta(U, shape1 = v_uc * (1 - p_nt) / (1 - p_at), shape2 = v_uc * p_nt / (1 - p_at)))
  D1 <- ifelse(D0 == 1, 1L, D1)
  D <- ((1 - Trt) * D0 + Trt * D1)[cluster_id]
  table(D0, D1)
  
  Y0 <- r_uy * qnorm(U, mean = 0, sd = sqrt(icc))[cluster_id] + 
          rnorm(clusters, mean = 0, sd = sqrt((1 - r_uy^2) * icc))[cluster_id] + 
            rnorm(sum(n_j), mean = 0, sd = sqrt(1 - icc))
  Y <- Y0 + delta * D  
  
  data.frame(cluster = cluster_id, Trt = as.integer(Trt[cluster_id]), D = D, Y = Y)
}

# dat <- r_cluster_trial(clusters, p_trt, delta, size_mean, p_nt, p_at, v_uc, r_uy, icc)

#------------------------------------------------------
# Model-fitting function
#------------------------------------------------------

iv_est <- function(dat) {
  require(AER)
  require(clubSandwich)
  lm_fit <- lm(Y ~ D, data = dat)
  lm_CR2 <- coef_test(lm_fit, cluster = dat$cluster, vcov = "CR2", test = c("z", "Satterthwaite"))["D",,drop=FALSE]
  rownames(lm_CR2) <- NULL
  
  iv_fit <- ivreg(Y ~ D | Trt, data = dat)
  
  iv_CR0 <- coef_test(iv_fit, cluster = dat$cluster, vcov = "CR0", test = c("z", "Satterthwaite"))["D",,drop=FALSE]
  rownames(iv_CR0) <- NULL
    
  iv_CR2 <- coef_test(iv_fit, cluster = dat$cluster, vcov = "CR2", test = c("z", "Satterthwaite"))["D",,drop=FALSE]
  rownames(iv_CR2) <- NULL  

  data.frame(type = c("lm-CR2","iv-CR0","iv-CR2"), rbind(lm_CR2, iv_CR0, iv_CR2))
}

# iv_est(dat)

#------------------------------------------------------
# Calculate performance measures
#------------------------------------------------------

check_pvals <- function(x, alpha) {
  res <- lapply(alpha, function(a) mean(x < a))
  res <- as.data.frame(res)
  names(res) <- paste("alpha", alpha, sep = "_")
  list(res)
}

# check_pvals(x = runif(10000), alpha = c(0.01, 0.05, 0.10))

#-----------------------------------------------------------
# Simulation Driver
#-----------------------------------------------------------

simulate_IV <- function(replicates, clusters, p_trt = 0.5, size_mean = 5, 
                        p_nt = 0.1, p_at = p_nt, v_uc = 2, r_uy = 0.8, icc = 0.2,
                        alpha = c(.01, .05), seed = NULL) {
  require(dplyr) 
  require(tidyr)
  if (!is.null(seed)) set.seed(seed)
  
  reps <- replicate(replicates, {
    dat <- r_cluster_trial(clusters = clusters, p_trt = p_trt, size_mean = size_mean, 
                           p_nt = p_nt, p_at = p_at, v_uc = v_uc,
                           r_uy = r_uy, icc = icc)
    iv_est(dat)
  }, simplify = FALSE)
  
  bind_rows(reps) %>%
    group_by(type) %>%
    mutate(
      V = SE^2,
      p_z = ifelse(is.na(p_z), 1, p_z),
      p_Satt = ifelse(is.na(p_Satt), 1, p_Satt)
    ) %>%
    summarise(
      na_pct = mean(is.na(beta)),
      E_beta = mean(beta, na.rm = TRUE),
      V_beta = var(beta, na.rm = TRUE),
      E_var = mean(V, na.rm = TRUE),
      V_var = var(V, na.rm = TRUE),
      z = check_pvals(p_z, alpha = alpha),
      Satt = check_pvals(p_Satt, alpha = alpha)
    ) %>%
    gather("test","reject", z, Satt) %>%
    unnest(reject)
}

# simulate_IV(replicates = 10, clusters = 20, p_trt = 0.3, 
#             p_nt = 0.25, p_at = 0.25, v_uc = 5, seed = 749088071)


#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20171128)

design_factors <- list(
  clusters = seq(20, 60, 10), 
  p_trt = c(.3, .5), 
  icc = 0.2, 
  v_uc = c(1, 2, 5),
  p_nt = seq(0, 0.25, 0.05)
)

params <- expand.grid(design_factors)
params$replicates <- 5000
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

warnings()

parallel::stopCluster(clust)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, file = "auxilliary/Clustered IV Simulation Results.Rdata")


