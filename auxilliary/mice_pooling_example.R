library(dplyr)
library(clubSandwich)
library(mice)
library(mlmRev)
rm(list=ls())

data(bdf)

# create missing values to impute

bdf$langPRET[sample(1:nrow(bdf),75)] = NA
bdf_subset = bdf[,c("schoolNR", "IQ.verb", "IQ.perf", "sex", "Minority", "langPRET", "langPOST")]


# impute missing values

Impute_bdf <- mice(bdf_subset, m=28, meth="norm.nob", seed=24)


# fit and pool results based on homoskedastic variance-covariance estimates

models_fragile <- with(data = Impute_bdf, lm(langPOST ~ langPRET + Minority + sex + IQ.perf + IQ.verb))

models_fragile %>%
  pool() %>%
  summary()


# fit results with clubSandwich standard errors

models_robust <- with(data=Impute_bdf, 
                      lm(langPOST ~ langPRET + Minority + sex + IQ.perf + IQ.verb) %>% 
                        coef_test(cluster=bdf_subset$schoolNR, vcov="CR2", test="Satterthwaite")
                      ) 


# pool results with clubSandwich standard errors

models_robust$analyses %>%
  
  # add coefficient names as a column
  lapply(function(x) {
    x$coef <- row.names(x)
    x
  }) %>%
  bind_rows() %>%
  as.data.frame() %>%
  
  # summarize by coefficient
  group_by(coef) %>%
  summarise(
    m = n(),
    V_betw = var(beta),
    est = mean(beta),
    V_com = mean(SE^2),
    df_com = mean(df)
  ) %>%
  
  mutate(
    
    # calculate intermediate quantities to get df
    V_total = V_com + (1 + 1 / m) * V_betw,
    fmi = (1 + 1 / m) * V_betw / V_total,
    df_m = (m - 1) / fmi^2,
    lambda = (df_com + 1) / (df_com + 3),
    df_obs = lambda * df_com * (1 - fmi),
    df = 1 / (1 / df_m + 1 / df_obs),
    
    # calculate summary quantities for output
    se = sqrt(V_total),
    t = est / se,
    p_val = 2 * pt(abs(t), df = df, lower.tail = FALSE),
    crit = qt(0.975, df = df),
    lo95 = est - se * crit,
    hi95 = est + se * crit
  ) %>%
  select(coef, est, t, df_com, df, p_val, lo95, hi95, fmi)

# note I left df_com (average of the complete data degrees of freedom) in there for comparison purposes