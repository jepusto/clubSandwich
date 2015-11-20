library(tidyr)
library(dplyr)
library(ggplot2)

rm(list=ls())
load("paper_ClusterRobustTesting/R/Panel simulation results.Rdata")


results <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
  UB <- alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)
})
head(results)

results_long <- gather(results, "hypothesis", "reject", t_B:F_all)  

filter(results_long, test %in% c("CR1 Naive-F", "CR2 HTZ")) %>%
  mutate(test = ifelse(test=="CR1 Naive-F", "Naive", "AHZ")) %>%
  spread(test, reject) ->
  results_compare


# ICC = 0

filter(results_compare, alpha==0.05 & icc==0) %>%
  group_by(design, m, n, hypothesis, alpha) %>%
  summarise(AHZ = mean(AHZ), Naive = mean(Naive), iterations = sum(iterations)) %>%
  mutate(UB = alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)) %>%
ggplot(aes(Naive, AHZ, color = factor(design))) + 
  geom_point() + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  geom_abline(slope = 1) + 
  facet_wrap(~ hypothesis, scales = "free") + 
  labs(color = "design") + 
  theme_minimal()

# ICC > 0

filter(results_compare, alpha==0.05 & icc > 0) %>%
  ggplot(aes(Naive, AHZ, color = factor(design))) + 
  geom_point() + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  geom_abline(slope = 1) + 
  facet_grid(hypothesis ~ trt_var, scales = "free", labeller = "label_both") + 
  labs(color = "design") + 
  theme_minimal()

