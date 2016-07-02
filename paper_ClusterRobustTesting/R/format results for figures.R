load("R/Panel simulation results.Rdata")

results <- within(results, {
  constraints <- NULL
  seed <- NULL
  design <- names(design)
})
head(results)

gather(results, "alpha", "reject", alpha0.005, alpha0.01, alpha0.05, alpha0.1) %>%
  mutate(alpha = as.numeric(substring(alpha, 6))) ->
  results_long_all

filter(results_long_all, test %in% c("CR1 Naive-F", "CR2 HTZ") 
       & design %in% c("CR-balanced","CR-unbalanced","DD-balanced-unbalanced",
                       "DD-unbalanced-unbalanced","RB-balanced-unequal","RB-unbalanced-unequal")) %>%
  within({
    test <- factor(ifelse(test=="CR1 Naive-F", "Standard", "AHT"), levels = c("Standard","AHT"))
    UB <- alpha + qnorm(0.95) * sqrt(alpha * (1 - alpha) / iterations)
    m_fac <- paste0("m = ",m)
    q <- hypothesis
    levels(q) <- c(1,1,2,3,3,6)
    design <- str_to_upper(str_extract(design, "[A-Z]+-[a-z]"))
    q_alpha <- paste0("q = ", q, ", alpha = ", alpha)
    alpha_q <- paste0("alpha = ", alpha, ", q = ", q)
    alpha_m <- paste0("alpha = ", alpha, ", m = ", m)
    test_q <- paste0(test, " test, q = ", q)
    m_test <- factor(paste0(test, " test, m = ", m))
  }) ->
  results_long
head(results_long)

select(results_long, -df) %>%
  spread(test, reject) ->
  results_compare
head(results_compare)


filter(results_long, n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-B") %>%
  mutate(reject=0) -> 
  zeros_long

filter(results_compare, m==15 & n==18 & icc==0.05 & 
         trt_var==0 & rho==0.2 & design=="CR-B") %>%
  mutate(Naive = 0, AHT = 0) -> 
  zeros_compare


alpha_val <- 0.4

breaks_cut <- function(alpha) {
  function(limits) { c(pretty(limits, 4),alpha)}
}


# Figure 1

filter(results_long, alpha == .05) %>%
  ggplot(aes(m_fac, reject)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val, fill = "blue", color = "blue") + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ test_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate") + 
  theme_bw() + theme(legend.position = "none")

# Figure 2

filter(results_long, alpha==0.05 & m==30) %>%
  ggplot(aes(design, reject, color = factor(design), fill = factor(design))) + 
  geom_boxplot(coef = Inf, alpha = alpha_val) + 
  geom_blank(data = zeros_long) + 
  geom_hline(aes(yintercept = alpha)) + 
  geom_hline(aes(yintercept = UB), linetype = "dashed") + 
  facet_wrap(~ test_q, ncol = 4, scales = "free") + 
  labs(x = NULL, y = "Rejection rate", color = "q", fill = "q") + 
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1))

# Figure 3

filter(results_long, alpha==0.05 & test=="AHT" & icc == 0.25 & rho == 0.8 & trt_var == 0.04) %>%
  mutate(m_lab = paste("m =", m)) %>%
  ggplot(aes(design, df, fill = design)) + 
  geom_boxplot(coef = Inf, alpha = alpha_val) + 
  facet_wrap(~ m_lab, scales = "free") + 
  labs(x = "Study design", y = "Denominator degrees of freedom") + 
  theme_bw() + theme(legend.position = "none")

