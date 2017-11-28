library(tidyverse)

#--------------------------------------------------------
# Analyze clustered simulation results
#--------------------------------------------------------
rm(list=ls())
load("auxilliary/Clustered IV Simulation Results.Rdata")

results$compliance <- with(results, 1 - 2 * p_nt)

results %>%
  filter(is.na(alpha_0.05)) %>%
  View()

results %>%
  filter(type == "lm-CR2", test == "Satt") %>%
  ggplot(aes(compliance, alpha_0.05, color = factor(v_uc), shape =  factor(v_uc), linetype =  factor(v_uc))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Number of clusters", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")

results %>%
  filter(type != "lm-CR2") %>%
  ggplot(aes(compliance, alpha_0.05, color = interaction(type,test), group = interaction(type, test, factor(v_uc)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Number of clusters", y = "Rejection rate at alpha = .05", color = "test") + 
  theme(legend.position = "bottom")

results %>%
  filter(type != "lm-CR2") %>%
  ggplot(aes(compliance, alpha_0.01, color = interaction(type,test), group = interaction(type, test, factor(v_uc)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .01, linetype = "dashed") + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Number of clusters", y = "Rejection rate at alpha = .05", color = "test") + 
  theme(legend.position = "bottom")

#--------------------------------------------------------
# Analyze multi-site simulation results
#--------------------------------------------------------

rm(list=ls())
load("auxilliary/IV Multisite Simulation Results.Rdata")

results %>%
  filter(is.na(alpha_0.05))


