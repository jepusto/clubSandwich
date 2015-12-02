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
  within({
    test <- ifelse(test=="CR1 Naive-F", "Naive", "AHZ")
    UB <- alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)
    m_fac <- paste0("m = ",m)
    q <- hypothesis
    levels(q) <- c(1,1,2,3,3,6)
    q_alpha <- paste0("q = ", q, ", alpha = ", alpha)
    alpha_q <- paste0("alpha = ", alpha, ", q = ", q)
    alpha_m <- paste0("alpha = ", alpha, ", m = ", m)
  }) ->
  results_long
head(results_long)

select(results_long, -df) %>%
  spread(test, reject) ->
  results_compare
head(results_compare)


filter(results_long, n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-balanced") %>%
  mutate(reject=0) -> 
  zeros_long

filter(results_compare, m==15 & n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-balanced") %>%
  mutate(Naive = 0, AHZ = 0) -> 
  zeros_compare

#--------------------------
# Overview
#--------------------------

# Take-away: Naive-F = bad, AHZ = rad!
# Take-away: Naive-F performance degrades (becomes liberal) for larger q
# Take-away: AHZ test becomes more conservative

# boxplots

ggplot(results_long, aes(m_fac, reject, fill = test)) + 
  geom_boxplot() + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ alpha_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate", fill = "Test") + 
  theme_bw() + theme(legend.position = "bottom")

# head-to-head, separate plots by alpha and q

ggplot(results_compare, aes(Naive, AHZ, color = m_fac, shape = m_fac)) + 
  geom_point() + 
  geom_blank(data = zeros_compare) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  facet_wrap(~ alpha_q, ncol = 4, scales = "free") + 
  labs(color = NULL, shape = NULL) + 
  theme_bw()

# head-to-head, separate plots by alpha and m

ggplot(results_compare, aes(Naive, AHZ, color = factor(q), shape = factor(q))) + 
  geom_point() + 
  geom_blank(data = zeros_compare) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  facet_wrap(~ alpha_m, ncol = 3, scales = "free") + 
  labs(color = NULL, shape = NULL) + 
  theme_bw()

#--------------------------
# Balance/leverage
#--------------------------

# Take-away: Actual rejection rates are strongly affected by balance

# alpha = .05, dotplot
filter(results_long, alpha==0.05) %>%
  ggplot(aes(design, reject, color = factor(q), shape = factor(q))) + 
  geom_point() + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ m, labeller = "label_both", scales = "free") + 
  labs(x = NULL, y = "Rejection rate", color = "q", shape = "q") + 
  theme_bw() + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# alpha = .05, boxplot excluding q = 6

filter(results_long, alpha==0.05 & q != 6) %>%
ggplot(aes(design, reject, fill = design)) + 
  geom_boxplot() + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ m, labeller = "label_both", scales = "free") + 
  labs(x = NULL, y = "Rejection rate", fill = NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

# alpha = .05, boxplot for q = 6

filter(results_long, alpha==0.05 & q == 6) %>%
  ggplot(aes(design, reject, fill = design)) + 
  geom_boxplot() + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ m, labeller = "label_both", scales = "free") + 
  labs(x = NULL, y = "Rejection rate", fill = NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

#--------------------------
# Model mis-specification
#--------------------------

# Take-away: Model mis-specification doesn't matter for AHZ. 

filter(results_long, alpha==0.05 & test=="AHZ") %>%
  ggplot(aes(factor(trt_var), reject, fill = factor(icc))) + 
  geom_boxplot() + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(m ~ q, labeller = "label_both", scales = "free") + 
  labs(x = "treatment effect variance", y = "Rejection rate", fill = "intra-class correlation") + 
  theme_bw() + theme(legend.position = "bottom")

