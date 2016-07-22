library(foreign)
library(dplyr)

MortalityRates <- 
  read.dta("http://masteringmetrics.com/wp-content/uploads/2015/01/deaths.dta", 
           convert.factors=FALSE) %>%
  filter(dtype %in% c(1,2,3,6) & agegr==2) %>%
  mutate(cause = factor(dtype, labels = c("All","Motor Vehicle","Suicide","Internal"))) %>%
  select(-dtype, -agegr, -age, -legal1820)

save(MortalityRates, file = "data/MortalityRates.RData")
dim(MortalityRates)
head(MortalityRates)