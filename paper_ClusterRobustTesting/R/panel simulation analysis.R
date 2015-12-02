library(tidyr)
library(dplyr)
library(ggplot2)

rm(list=ls())
load("paper_ClusterRobustTesting/R/Panel simulation results.Rdata")


results <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
})
head(results)

gather(results, "alpha", "reject", alpha1:alpha10) %>%
  mutate(alpha = as.numeric(substring(alpha, 6)) / 100) ->
  results_long_all

filter(results_long_all, test %in% c("CR1 Naive-F", "CR2 HTZ")) %>%
  mutate(test = ifelse(test=="CR1 Naive-F", "Naive", "AHZ"),
         UB = alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations),
         q = switch(hypothesis, t_B = 1, t_C = 1, F_1 = 2, F_B = 3, F_C = 3, F_all = 6),
         q_alpha = paste0("q = ", q, "alpha = ", alpha)) ->
  results_long

filter(results_long_all, test %in% c("CR1 Naive-F", "CR2 HTZ")) %>%
  mutate(test = ifelse(test=="CR1 Naive-F", "Naive", "AHZ"),
         UB = alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)) %>%
  select(-df) %>%
  spread(test, reject) ->
  results_compare

#--------------------------
# Overview
#--------------------------

# Take-away: Naive-F = bad, AHZ = rad!
filter(results_long, n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-balanced") %>%
  select(m, reject, test, hypothesis) %>%
  mutate(reject=0) -> 
  zeros

ggplot(results_long, aes(factor(m), reject, fill = test)) + 
  geom_boxplot() + 
  geom_blank(data = zeros) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ hypothesis_alpha, nrow = 3, scales = "free") + 
  labs(x = "Clusters (m)", y = "Rejection rate", fill = "Test") + 
  theme_bw() + theme(legend.position = "bottom")

  
#--------------------------
# Balance/leverage
#--------------------------

# Take-away: Actual rejection rates are strongly affected by balance

#--------------------------
# Number of constraints
#--------------------------

# Take-away: Naive-F performance degrades (becomes liberal) for larger q
# AHZ test becomes more conservative

#--------------------------
# Model mis-specification
#--------------------------

# Take-away: Model mis-specification doesn't matter very much for AHZ. 

#--------------------------
# Head-to-head comparisons
#--------------------------

# alpha = .01

filter(results_compare, alpha==0.01) %>%
  ggplot(aes(Naive, AHZ, color = factor(design))) + 
  geom_point() + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  geom_abline(slope = 1) + 
  facet_grid(m ~ hypothesis, scales = "free") + 
  labs(color = "design") + 
  theme_minimal()


# alpha = .05

filter(results_compare, alpha==0.05) %>%
  ggplot(aes(Naive, AHZ, color = factor(design))) + 
  geom_point() + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  geom_abline(slope = 1) + 
  facet_grid(m ~ hypothesis, scales = "free") + 
  labs(color = "design") + 
  theme_minimal()


# alpha = .10

filter(results_compare, alpha==0.10) %>%
  ggplot(aes(Naive, AHZ, color = factor(design))) + 
  geom_point() + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  geom_abline(slope = 1) + 
  facet_grid(m ~ hypothesis, scales = "free") + 
  labs(color = "design") + 
  theme_minimal()
