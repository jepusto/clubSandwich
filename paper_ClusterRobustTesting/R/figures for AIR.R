rm(list=ls())
setwd("paper_ClusterRobustTesting")
source("R/format results for figures.R")

# standard test

filter(results_long, alpha == .05 & test == "CR1 standard") %>%
  mutate(q_fac = paste("q =",q,ifelse(q==1, "(t-test)", "(F-test)"))) %>%
  ggplot(aes(m_fac, reject)) + 
  geom_boxplot(coef = Inf,  fill = "blue") + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ q_fac, ncol = 2, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")
ggsave("AIR standard test.wmf", width = 8, height = 5)

# AHT test
filter(results_long, alpha == .05 & test %in% c("CR1 standard", "CR2 AHT")) %>%
  ggplot(aes(m_fac, reject, fill = test_lab)) + 
  geom_boxplot(coef = Inf) + 
  geom_blank(data = filter(zeros_long, test %in% c("CR1 standard", "CR2 AHT"))) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  scale_y_continuous(breaks = breaks_cut(.05)) + 
  facet_wrap(~ test_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")
ggsave("AIR AHT test.wmf", width = 10, height = 5)
