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
  results_long

filter(results_long, test %in% c("CR1 Naive-F", "CR2 HTZ")) %>%
  mutate(test = ifelse(test=="CR1 Naive-F", "Naive", "AHZ"),
         UB = alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)) %>%
  select(-df) %>%
  spread(test, reject) ->
  results_compare


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
