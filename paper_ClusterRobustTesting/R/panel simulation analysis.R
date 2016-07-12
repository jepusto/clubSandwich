library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

rm(list=ls())
load("paper_ClusterRobustTesting/R/Panel simulation results.Rdata")

results <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
})
head(results)

gather(results, "alpha", "reject", alpha0.005, alpha0.01, alpha0.05, alpha0.1) %>%
  mutate(alpha = as.numeric(substring(alpha, 6))) ->
  results_long_all

filter(results_long_all, test %in% c("CR1 Naive-F", "CR1 HTZ","CR2 Naive-F","CR2 HTZ", "CR2A HTZ")) %>%
  within({
    test <- gsub("Naive-F","standard", gsub("HTZ","AHT", test))
    test_lab <- ifelse(str_detect(test, "standard"), "Standard", ifelse(str_detect(test, "AHT"), "AHT",test))
    UB <- alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)
    m_fac <- paste0("m = ",m)
    q <- hypothesis
    levels(q) <- c(1,1,2,3,3,6)
    design <- str_to_upper(str_extract(design, "[A-Z]+-[a-z]"))
    q_alpha <- paste0("q = ", q, ", alpha = ", alpha)
    alpha_q <- paste0("alpha = ", alpha, ", q = ", q)
    alpha_m <- paste0("alpha = ", alpha, ", m = ", m)
    test_q <- paste0(test_lab, " test, q = ", q)
    m_test <- factor(paste0(test, " test, m = ", m))
  }) ->
  results_long
head(results_long)

select(results_long, -df, -m_test) %>%
  mutate(test = gsub(" ","",test)) %>%
  spread(test, reject) ->
  results_compare
head(results_compare)


filter(results_long, n==18 & icc==0.05 & trt_var==0 & rho==0.2 & design=="CR-B") %>%
  mutate(reject=0) -> 
  zeros_long

filter(results_compare, m==15 & n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-balanced") %>%
  mutate(CR1adhoc = 0, CR1AHT = 0, CR2AAHT = 0, CR2adhoc = 0, CR2AHT = 0) -> 
  zeros_compare

alpha_val <- 0.4

breaks_cut <- function(alpha) {
  function(limits) { 
    if (max(limits) < 6 * alpha) {
      c(pretty(limits, 4), alpha) 
    } else {
      pretty(limits, 4)
    }
  }
}

#--------------------------
# Overview
#--------------------------

# Take-away: Naive-F = bad, AHZ = rad!
# Take-away: Naive-F performance degrades (becomes liberal) for larger q
# Take-away: AHZ test becomes more conservative

# boxplots
filter(results_long, alpha == .05 & test %in% c("CR1 standard", "CR2 AHT")) %>%
  ggplot(aes(m_fac, reject)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val, fill = "grey") + 
  geom_blank(data = filter(zeros_long, test %in% c("CR1 standard", "CR2 AHT"))) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  scale_y_continuous(breaks = breaks_cut(.05)) + 
  facet_wrap(~ test_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")

# head-to-head, separate plots by alpha and q

ggplot(results_compare, aes(CR1AHT, CR2AHT, color = m_fac, shape = m_fac)) + 
  geom_point() + 
  geom_blank(data = zeros_compare) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  facet_wrap(~ alpha_q, ncol = 4, scales = "free") + 
  labs(color = NULL, shape = NULL) + 
  theme_bw()


ggplot(results_compare, aes(CR2AAHT, CR2AHT, color = m_fac, shape = m_fac)) + 
  geom_point() + 
  geom_blank(data = zeros_compare) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_vline(aes(xintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_vline(aes(xintercept = UB), linetype = "dashed") + 
  facet_wrap(~ alpha_q, ncol = 4, scales = "free") + 
  labs(color = NULL, shape = NULL) + 
  theme_bw()

#--------------------------
# Balance/leverage
#--------------------------

# Take-away: Actual rejection rates are strongly affected by balance

# alpha = .05, dotplot

filter(results_long, alpha==0.05 & m==30 & test %in% c("CR1 standard", "CR2 AHT")) %>%
  ggplot(aes(design, reject)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val, fill = "grey") + 
  geom_blank(data = filter(zeros_long, test %in% c("CR1 standard", "CR2 AHT"))) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  scale_y_continuous(breaks = breaks_cut(.05)) + 
  facet_wrap(~ test_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

# alpha = .05, squished boxplots

# alpha = .05, boxplot excluding q = 6

filter(results_long, alpha==0.05 & q != 6) %>%
ggplot(aes(design, reject, fill = design)) + 
  geom_boxplot(coef = Inf) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ m, labeller = "label_both", scales = "free") + 
  labs(x = NULL, y = "Rejection rate", fill = NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

# alpha = .05, boxplot for q = 6

filter(results_long, alpha==0.05 & q == 6) %>%
  ggplot(aes(design, reject, fill = design)) + 
  geom_boxplot(coef = Inf) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(test ~ m, labeller = "label_both", scales = "free") + 
  labs(x = NULL, y = "Rejection rate", fill = NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

filter(results_long, alpha==.05 & q==6 & m==15 & test %in% c("CR2 AHT","CR2A AHT")) %>%
  select(design, n, icc, trt_var, rho, ar, test, df, reject)

#--------------------------
# Model mis-specification
#--------------------------

# Take-away: Model mis-specification doesn't matter for AHZ. 

filter(results_long, alpha==0.05 & test=="CR2 AHT") %>%
  ggplot(aes(factor(trt_var), reject, fill = factor(icc), color = factor(icc))) + 
  geom_boxplot(coef = Inf, alpha = alpha_val) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_grid(m ~ q, labeller = "label_both", scales = "free") + 
  labs(x = "treatment effect variance", y = "Rejection rate", fill = "intra-class correlation", color = "intra-class correlation") + 
  theme_bw() + theme(legend.position = "bottom")

# Look at whether adding icc/correlation structure improves rejection rate accuracy for very small sample sizes

#--------------------------
# Degrees of freedom
#--------------------------

filter(results_long, alpha==0.05 & test=="CR2 AHT" & icc == 0.25 & rho == 0.8 & trt_var == 0.04) %>%
  mutate(m_lab = paste("m =", m)) %>%
  ggplot(aes(design, df)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val, fill = "grey") + 
  facet_wrap(~ m_lab, scales = "free") + 
  labs(x = "Study design", y = "Denominator degrees of freedom") + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

filter(results_long, alpha==0.05 & test %in% c("CR1 AHT", "CR2 AHT", "CR2A AHT") & icc == 0.25 & rho == 0.8 & trt_var == 0.04) %>%
  mutate(m_lab = paste("m =", m)) %>%
  ggplot(aes(design, df, fill = test, color = test)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val) + 
  facet_wrap(~ m_lab, scales = "free") + 
  labs(x = "Study design", y = "Denominator degrees of freedom") + 
  theme_bw() + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))




#--------------------------
# CR2A figures
#--------------------------

filter(results_long, alpha==0.01 & 
         test %in% c("CR2 AHT","CR2A AHT") &
         design %in% c("CR-balanced","CR-unbalanced","DD-balanced","DD-unbalanced")) %>%
  mutate(test = ifelse(test=="CR2 AHT", "Accounting for absorption", "Ignoring absorption")) %>%
  group_by(m_fac, design, test, q, alpha, UB) %>%
  summarise(min = min(reject), lower = quantile(reject, .25), middle = median(reject), 
            upper = quantile(reject, .75), max = max(reject)) %>%
  ggplot(aes(m_fac, fill = test, color = test)) + 
  geom_boxplot(aes(ymin = min, lower = lower, middle = middle, upper = upper, ymax = max), 
               stat = "identity", alpha = alpha_val) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  scale_y_continuous(breaks = breaks_cut(.05)) + 
  scale_fill_brewer(type = "qual", palette = 6) + 
  facet_grid(q ~ design, scales = "free_y") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "bottom")


zeros_long_temp <- 
  filter(zeros_long, test %in% c("CR2 AHT","CR2A AHT")) %>%
  mutate(test = ifelse(test=="CR2 AHT", "Accounting for absorption", "Ignoring absorption"))

filter(results_long, 
         test %in% c("CR2 AHT","CR2A AHT") &
         design %in% c("CR-balanced","CR-unbalanced","DD-balanced","DD-unbalanced")) %>%
  mutate(test = ifelse(test=="CR2 AHT", "Accounting for absorption", "Ignoring absorption")) %>%
  group_by(m_fac, test, q, alpha, q_alpha, UB) %>%
  summarise(min = min(reject), lower = quantile(reject, .25), middle = median(reject), 
            upper = quantile(reject, .75), max = max(reject)) %>%
  ggplot(aes(m_fac, fill = test, color = test)) + 
  geom_boxplot(aes(ymin = min, lower = lower, middle = middle, upper = upper, ymax = max), 
               stat = "identity", alpha = alpha_val) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  geom_blank(data = zeros_long_temp) + 
  scale_y_continuous(breaks = breaks_cut(.05)) + 
  scale_fill_brewer(type = "qual", palette = 6) + 
  facet_grid(alpha ~ q, scales = "free_y", labeller = "label_both") + 
  labs(x = NULL, y = "Rejection rate", fill = "CR2 adjustment", color = "CR2 adjustment") + 
  theme_bw() + theme(legend.position = "bottom")
