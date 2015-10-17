setwd("paper_ClusterRobustTesting")
library(tidyr)
library(dplyr)
library(readstata13)
library(plm)
read.dta13("data/CPS_all_micro.dta", convert.factors=FALSE) %>%
  droplevels() ->
  all
names(all)
nrow(all) # 1925329 
nlevels(factor(all$statefip)) # 51 states
nlevels(factor(all$year)) # 36 years

all$female <- all$sex - 1

select(all, statefip, year, lnwage, yrseduc, age, age2, female) %>%
group_by(statefip, year) %>%
  summarise_each(funs(mean)) %>%
  mutate(year = factor(year)) ->
  CPS_agg_panel

plm_fit <- plm(lnwage ~ yrseduc + age + age2 + female, 
               data = CPS_agg_panel, index = c("statefip","year"))
summary(plm_fit)

# all$age <- as.numeric(paste(all$age)) + 1
# 
# panel2 <- aggregate(all,by = list(all$year,all$statefip),FUN = "mean",na.action=na.omit)
# panel3 <- panel2[c("year","Group.2","wtsupp","yrseduc","wage_per_hour","lnwage","age","age2","female")]
# 
# colnames(panel3) <- c("year","statefip","wtsupp","yrseduc","wage_per_hour","lnwage","age","age2","female")
# 
# panel3$age_s <- ave(panel3$age,panel3$statefip,FUN = mean)
# panel3$age2_s <- ave(panel3$age2,panel3$statefip,FUN = mean)
# panel3$female_s <- ave(panel3$female,panel3$statefip,FUN = mean)
# panel3$yrseduc_s <- ave(panel3$yrseduc,panel3$statefip,FUN = mean)
# panel3$lnwage_s <- ave(panel3$yrseduc,panel3$statefip,FUN= mean)
# 
# panel3$age_f <- panel3$age - panel3$age_s
# panel3$age2_f <- panel3$age2 - panel3$age2_s
# panel3$female_f <- panel3$female - panel3$female_s
# panel3$yrseduc_f <- panel3$yrseduc - panel3$yrseduc_s
# panel3$lnwage_f <- panel3$lnwage - panel3$lnwage_s

#####################
#Hausman Test Analysis
#####################
#FE model
FE <- lm(lnwage_f~  yrseduc_f + age_f + age2_f + female_f, data = panel3)

#Between model
panel4 <- aggregate(panel3, by = list(panel3$statefip), FUN = "mean", na.action = na.omit)
Bet <- lm(lnwage~ yrseduc + age + age2 + female, data = panel4 )

#Theta calculation (following definitions; dift than thru PLM?)
n <- 51 #number clusters
time <- 36 #number of time periods
k <- 4 #number of covariates
SigE2 <- sum((FE$resid)^2)/(n*(time - 1) - k) #within-var
SigB2 <- sum((Bet$resid)^2)/(n - k - 1) #between-var
SigU2 <- SigB2 - SigE2/time
theta <- 1 - sqrt(SigE2/(SigE2 + time*SigU2)) 
#theta <- 1 - sqrt(SigE2/(time*SigB2)) #equiv
theta #.4582

#Quasi-demeaning
panel3$constant_r <- (1 - theta)
panel3$age_r <- panel3$age - theta*panel3$age_s
panel3$age2_r <- panel3$age2 - theta*panel3$age2_s
panel3$female_r <- panel3$female - theta*panel3$female_s
panel3$yrseduc_r <- panel3$yrseduc - theta*panel3$yrseduc_s
panel3$lnwage_r <- panel3$lnwage - theta*panel3$lnwage_s

#RE model
RE <- lm(lnwage_r ~ -1 + constant_r + yrseduc_r + age_r + age2_r + female_r, 
         data = panel3)

#Testing FE vs RE: Alternative Hausman test
REfull <- lm(lnwage_r ~ yrseduc_r + age_r + age2_r + female_r 
             + yrseduc_s + age_s + age2_s + female_s, 
             data = panel3)

coef_test(REfull, vcov = "CR2", cluster = panel3$statefip)
coef_test(FE,vcov="CR2",cluster=panel3$statefip)
vcov_h <- vcovCR(REfull, cluster = panel3$statefip, type = "CR2", inverse_var = TRUE)
Wald_test(REfull, 6:9 , vcov = vcov_h, test = "All") 
