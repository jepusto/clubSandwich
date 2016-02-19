library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
rm(list=ls())
setwd("paper_ClusterRobustTesting")
source("R/format results for figures.R")

# standard test

filter(results_long, alpha == .05 & test == "Standard") %>%
  mutate(q_fac = paste("q =",q,ifelse(q==1, "(t-test)", "(F-test)"))) %>%
  ggplot(aes(m_fac, reject)) + 
  geom_boxplot(coef = Inf,  fill = "blue") + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ q_fac, ncol = 2, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")
ggsave("standard test.wmf", width = 8, height = 5)

# BRL + Satterthwaite t-test

filter(results_long, alpha == .05 & q == 1) %>%
  mutate(test_lab = ifelse(test=="AHT","BRL + Satterthwaite t-test","Standard t-test")) %>%
  ggplot(aes(m_fac, reject, fill = test)) + 
  geom_boxplot(coef = Inf) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ test_lab) + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")
ggsave("BRL + Satt.wmf", width = 7, height = 3.5)

# AHT test

filter(results_long, alpha == .05 & q != 1) %>%
  mutate(q_fac = paste("q =",q, "(F-test)")) %>%
  droplevels() %>%
  ggplot(aes(m_fac, reject, fill = test)) + 
  geom_boxplot(coef = Inf) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ q_fac, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")
ggsave("AHT test.wmf", width = 8, height = 5)
