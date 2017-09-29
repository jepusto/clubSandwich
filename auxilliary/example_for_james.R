library(mice)
library(mlmRev)
rm(list=ls())
# source('auxilliary/mice_external.R')
data(bdf)

# create missing values to impute
bdf$langPRET[sample(1:nrow(bdf),75)] = NA
bdf_subset = bdf[,c("schoolNR", "IQ.verb", "IQ.perf", "sex", "Minority", "langPRET", "langPOST")]

Impute_bdf <- mice(bdf_subset, m=28, meth="norm.nob", seed=24)

models <- with(data=Impute_bdf, 
               coef_test(lm(langPOST ~ langPRET + Minority + sex + IQ.perf + IQ.verb), 
                         cluster=bdf_subset$schoolNR, vcov="CR2", 
                         test="Satterthwaite")) 

library(dplyr)

models$analyses %>%
  lapply(function(x) {
    x$coef <- row.names(x)
    x
  }) %>%
  bind_rows() %>%
  as.data.frame() %>%
  group_by(coef) %>%
  summarise(
    m = n(),
    V_betw = var(beta),
    beta = mean(beta),
    V_com = mean(SE^2),
    df_com = mean(df)
  ) %>%
  mutate(
    V_total = V_com + (1 + 1 / m) * V_betw,
    f = (1 + 1 / m) * V_betw / V_total,
    df_m = (m - 1) / f^2,
    lambda = (df_com + 1) / (df_com + 3),
    df_obs = lambda * df_com * (1 - f),
    df = 1 / (1 / df_m + 1 / df_obs)
  )
