library(tidyverse)

#--------------------------------------------------------
# Analyze clustered simulation results
#--------------------------------------------------------
rm(list=ls())
load("auxilliary/Clustered IV Simulation Results.Rdata")

results$compliance <- with(results, 1 - 2 * p_nt)

results %>%
  filter(is.na(alpha_0.05))

# bias of OLS estimates

results %>%
  filter(type == "lm-CR2", test == "Satt") %>%
  ggplot(aes(compliance, E_beta, color = factor(v_uc), shape =  factor(v_uc), linetype =  factor(v_uc))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Number of clusters", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")

# bias of IV estimates

results %>%
  filter(type == "iv-CR0", test == "Satt") %>%
  ggplot(aes(compliance, E_beta, color = factor(v_uc), shape =  factor(v_uc), linetype =  factor(v_uc))) + 
  geom_line() + geom_point() + 
  expand_limits(y = 0) + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")

# rejection rates of OLS tests

results %>%
  filter(type == "lm-CR2", test == "Satt") %>%
  ggplot(aes(compliance, alpha_0.05, color = factor(v_uc), shape =  factor(v_uc), linetype =  factor(v_uc))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")

# rejection rates of IV tests at alpha = .05

results %>%
  filter(type != "lm-CR2") %>%
  ggplot(aes(compliance, alpha_0.05, color = interaction(type,test), group = interaction(type, test, factor(v_uc)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .01", color = "test") + 
  theme(legend.position = "bottom")


# rejection rates of IV tests at alpha = .01

results %>%
  filter(type != "lm-CR2") %>%
  ggplot(aes(compliance, alpha_0.01, color = interaction(type,test), group = interaction(type, test, factor(v_uc)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .01, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(p_trt ~ clusters, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .01", color = "test") + 
  theme(legend.position = "bottom")

# For SREE 2019 abstract 

SREE_results <- 
  results %>%
  filter(
    (type == "iv-CR2" & test == "Satt") | (type == "iv-CR0" & test == "z")
  ) %>%
  mutate(test = recode(test, `Satt` = "Modified", `z` = "Conventional")) %>% 
  rename(`Allocation fraction` = p_trt)

ggplot(SREE_results, aes(compliance, alpha_0.05, color = test, linetype = factor(v_uc))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(`Allocation fraction` ~ clusters, labeller = "label_both") + 
  scale_color_brewer(type = "qual", palette = 6) + 
  scale_linetype(guide = "none") + 
  theme_light() +
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05", color = "CRVE estimator") + 
  theme(legend.position = "bottom", strip.text = element_text(color = "black"))

ggplot(SREE_results, aes(compliance, alpha_0.01, color = test, linetype = factor(v_uc))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .01, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(`Allocation fraction` ~ clusters, labeller = "label_both") + 
  scale_color_brewer(type = "qual", palette = 6) + 
  scale_linetype(guide = "none") + 
  theme_light() +
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .01", color = "CRVE estimator") + 
  theme(legend.position = "bottom", strip.text = element_text(color = "black"))


#--------------------------------------------------------
# Analyze multi-site simulation results
#--------------------------------------------------------

rm(list=ls())
load("auxilliary/Multisite IV Simulation Results.Rdata")

results <- 
  results %>%
  group_by(sites, delta_sd, compliance, comp_sd, type, test) %>%
  summarise_at(vars(E_beta, V_beta, E_var, V_var, alpha_0.01, alpha_0.05), .funs = mean)

# bias of OLS estimates

results %>%
  filter(type == "lm-CR2", test == "Satt") %>%
  ggplot(aes(compliance, E_beta, color = factor(comp_sd), shape =  factor(comp_sd), linetype =  factor(comp_sd))) + 
  geom_line() + geom_point() + 
  facet_grid(delta_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Bias of OLS estimates") + 
  theme(legend.position = "bottom")

# bias of IV estimates

results %>%
  filter(type %in% c("iv-one-CR0", "iv-many-CR0"), test == "Satt") %>%
  ggplot(aes(compliance, E_beta, color = type, group = interaction(type, factor(delta_sd)))) + 
  geom_line() + geom_point() + 
  expand_limits(y = 0) + 
  facet_grid(comp_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Bias of IV estimates", color = "test") + 
  theme(legend.position = "bottom")

# rejection rates of OLS tests

results %>%
  filter(type == "lm-CR2", test == "Satt") %>%
  ggplot(aes(compliance, alpha_0.05, color = factor(comp_sd), shape =  factor(comp_sd), linetype =  factor(comp_sd))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  facet_grid(delta_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05") + 
  theme(legend.position = "bottom")

# rejection rates of single-instrument tests at alpha = .05

results %>%
  filter(type %in% c("iv-one-CR0", "iv-one-CR2")) %>%
  ggplot(aes(compliance, alpha_0.05, color = interaction(type,test), group = interaction(type, test, factor(delta_sd)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(comp_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05", color = "test") + 
  theme(legend.position = "bottom")

# rejection rates of site-instrument tests at alpha = .05

results %>%
  filter(type %in% c("iv-many-CR0", "iv-many-CR2")) %>%
  ggplot(aes(compliance, alpha_0.05, color = interaction(type,test), group = interaction(type, test, factor(delta_sd)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .05, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(comp_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .05", color = "test") + 
  theme(legend.position = "bottom")


# rejection rates of single-instrument tests at alpha = .01

results %>%
  filter(type %in% c("iv-one-CR0", "iv-one-CR2")) %>%
  ggplot(aes(compliance, alpha_0.01, color = interaction(type,test), group = interaction(type, test, factor(delta_sd)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .01, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(comp_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .01", color = "test") + 
  theme(legend.position = "bottom")

# rejection rates of site-instrument tests at alpha = .01

results %>%
  filter(type %in% c("iv-many-CR0", "iv-many-CR2")) %>%
  ggplot(aes(compliance, alpha_0.01, color = interaction(type,test), group = interaction(type, test, factor(delta_sd)))) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept = .01, linetype = "dashed") + 
  expand_limits(y = 0) + 
  facet_grid(comp_sd ~ sites, labeller = "label_both") + 
  theme_light() + 
  labs(x = "Compliance rate", y = "Rejection rate at alpha = .01", color = "test") + 
  theme(legend.position = "bottom")
