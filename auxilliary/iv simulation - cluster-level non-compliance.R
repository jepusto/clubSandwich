rm(list = ls())

#------------------------------------------------------
# Data Generating Model
#------------------------------------------------------

# clusters <- 20
# trt_clusters <- 14
# delta <- 0
# size_mean <- 5
# c1 <- qnorm(0.1)
# c2 <- qnorm(0.8)
# c3 <- 0.7
# r_uy <- 0.8
# icc <- 0.3

r_cluster_trial <- function(clusters, trt_clusters, 
                            delta = 0, size_mean = 10, 
                            c1 = qnorm(0.1), c2 = qnorm(0.9), c3 = 0.7, 
                            r_uy = 0.8, icc = 0.3) {
  
  n_j <- 1 + rpois(clusters, lambda = size_mean)
  cluster_id <- rep(1:clusters, n_j)
  Trt <- (1:clusters) %in% (1:trt_clusters)
  
  U <- rnorm(clusters, mean = 0, sd = sqrt(icc))
  D0 <- rbinom(clusters, size = 1, prob = pnorm(c1 + c3 * U))
  D1 <- rbinom(clusters, size = 1, prob = pnorm(c2 + c3 * U))
  D1[D1 < D0] <- D0[D1 < D0]
  D <- ((1 - Trt) * D0 + Trt * D1)[cluster_id]
  
  Y0 <- rnorm(clusters, mean = r_uy * U, sd = sqrt((1 - r_uy^2) * icc))[cluster_id] + rnorm(sum(n_j), mean = 0, sd = sqrt(1 - icc))
  Y <- Y0 + delta * D  
  
  data.frame(cluster = cluster_id, Trt = as.integer(Trt[cluster_id]), D = D, Y = Y)
}

# dat <- r_cluster_trial(clusters, trt_clusters, delta, size_mean, c1, c2, c3, r_uy, icc)

#------------------------------------------------------
# Model-fitting function
#------------------------------------------------------

iv_est <- function(dat) {
  require(AER)
  require(clubSandwich)
  lm_fit <- lm(Y ~ D, data = dat)
  lm_CR2 <- coef_test(lm_fit, cluster = dat$cluster, vcov = "CR2", test = "Satterthwaite")["D",c("beta","SE","p_Satt"),drop=FALSE]
  rownames(lm_CR2) <- NULL
  names(lm_CR2)[3] <- "pval"

  iv_fit <- ivreg(Y ~ D | Trt, data = dat)
  iv_CR0 <- coef_test(iv_fit, cluster = dat$cluster, vcov = "CR0", test = "z")["D",,drop=FALSE]
  rownames(iv_CR0) <- NULL
  names(iv_CR0)[3] <- "pval"
  iv_CR2 <- coef_test(iv_fit, cluster = dat$cluster, vcov = "CR2", test = "Satterthwaite")["D",c("beta","SE","p_Satt"),drop=FALSE]
  rownames(iv_CR2) <- NULL
  names(iv_CR2)[3] <- "pval"
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

check_pvals(x = runif(10000), alpha = c(0.01, 0.05, 0.10))

#-----------------------------------------------------------
# Simulation Driver
#-----------------------------------------------------------

simulate_IV <- function(replicates, clusters, trt_clusters = clusters / 2, size_mean = 5, 
                        c1 = qnorm(0.1), c2 = -c1, c3 = 0.7, r_uy = 0.8, icc = 0.3,
                        alpha = c(.01, .05), seed = NULL) {
  require(dplyr) 
  require(tidyr)
  if (!is.null(seed)) set.seed(seed)
  
  reps <- replicate(replicates, {
    dat <- r_cluster_trial(clusters = clusters, trt_clusters = trt_clusters, 
                           size_mean = size_mean, c1 = c1, c2 = c2, c3 = c3,
                           r_uy = r_uy, icc = icc)
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
      reject = check_pvals(pval, alpha = alpha)
    ) %>%
    unnest(reject)
}

# simulate_IV(replicates = 100, clusters = 20, trt_clusters = 10, size_mean = 5,
#             c1 = qnorm(0.1), c2 = qnorm(0.9), c3 = 0.7,
#             r_uy = 0.8, icc = 0.3)


#-------------------------------------
# Experimental Design
#-------------------------------------
source_obj <- ls()

set.seed(20171119)

design_factors <- list(clusters = seq(20, 60, 10), c1 = qnorm(c(0.1, 0.2)), 
                       icc = c(0.1, 0.2, 0.3))
params <- expand.grid(design_factors)
params$replicates <- 50
params$seed <- round(runif(1) * 2^30) + 1:nrow(params)

# All look right?
lengths(design_factors)
nrow(params)
head(params)

#--------------------------------------------------------
# run simulations in parallel - mdply workflow
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

save(params, results, session_info, run_date, file = "auxilliary/IV Simulation Results.Rdata")


#--------------------------------------------------------
# Analyze simulation results
#--------------------------------------------------------
library(ggplot2)
rm(list=ls())
load("auxilliary/IV Simulation Results.Rdata")

results$compliance <- 1 - 2 * pnorm(results$c1)

ggplot(results, aes(clusters, alpha_0.05, color = type)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(compliance ~ icc, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Number of clusters", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")
