library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

load("R/Panel simulation results.Rdata")

results <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
})
head(results)

gather(results, "alpha", "reject", p0.005, p0.01, p0.05, p0.1) %>%
  mutate(alpha = as.numeric(substring(alpha, 2)),
         reject = ifelse(pNA==1, 0, reject)) ->
  results_long_all

filter(results_long_all, 
       test %in% c("CR1 Naive-F", "CR1 HTZ","CR2 Naive-F","CR2 HTZ", "CR2A HTZ") & 
         design %in% c("CR-balanced","CR-unbalanced","RB-balanced","RB-unbalanced","DD-balanced","DD-unbalanced")) %>%
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

# Figure 1

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

filter(results_long, m == 15 & test == "CR2 AHT") %>%
  group_by(alpha) %>%
  summarise(reject = max(reject)) ->
  AHT_max_reject

filter(results_long, alpha == .05 & test == "CR1 standard") %>%
  group_by(m, q) %>%
  summarise(reject = max(reject)) %>%
  spread(m, reject) ->
  standard_max_reject

# Figure 2

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

filter(results_long, m == 30 & test == "CR2 AHT" & alpha == .05) %>%
  summarise(min = min(reject), max = max(reject)) ->
  AHT_range_m30

# Figure 3

filter(results_long, alpha==0.05 & test=="CR2 AHT" & icc == 0.25 & rho == 0.8 & trt_var == 0.04) %>%
  mutate(m_lab = paste("m =", m)) %>%
  ggplot(aes(design, df)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val, fill = "grey") + 
  facet_wrap(~ m_lab, scales = "free") + 
  labs(x = "Study design", y = "Denominator degrees of freedom") + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

