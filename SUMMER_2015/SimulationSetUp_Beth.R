##simulation design ideas 

##-----------------------
#Basics
##-----------------------
#Always absorb cluster FE
#W = I
#m = 5, 10, 20, 40, 100
#n = 50, 100, 500 (bigger than MA sims)

#Types of tests
  #Hausman type tests (CPS balanced panel)
    #1-covariate w/in vs between
      #2 types: dummy/ continuous
    #Multi-covariates w/in vs between
      #Based on real data?

  #Treatment heterogeneity: Policy var x covariate
    #Policy dummy (assigned within-cluster)
      #Balanced vs not 
    #Covariates
      #within-covariate x policy
      #between-covariate x policy

##Keep in mind
  #From Tipton (2014):
    #cluster-level covariates have smallest DF
      #for dummies, clearly a function of prop 1's
      #for N(0,1), still closer to m/4
    #observation-level covariates 
      #holding m constant, as obs increases, df increases
      #nearly always > than df for cluster-level

##-----------------------------
#Explore CPS data 
##-----------------------------
getwd()
setwd("/Users/ltipton/Dropbox/Publications/Under Review/JEBS Tipton-Pustejovsky 2015/brlinr3/")
all <- read.dta13("data/HausmanTest/CPS_all_micro.dta")
sum(table(table(all$statefip))[-c(1)]) #51 states, with between 519 - 5866 observations
sum(table(table(all$year))) #36 years, with between 35189-75990 observations

all$female <- ifelse(all$sex == "Female",1,0)
all$age <- as.numeric(paste(all$age)) + 1

panel2 <- aggregate(all,by = list(all$year,all$statefip),FUN = "mean",na.action=na.omit)
panel3 <- panel2[c("year","Group.2","wtsupp","yrseduc","wage_per_hour","lnwage","age","age2","female")]

colnames(panel3) <- c("year","statefip","wtsupp","yrseduc","wage_per_hour","lnwage","age","age2","female")

panel3$age_s <- ave(panel3$age,panel3$statefip,FUN = mean)
panel3$age2_s <- ave(panel3$age2,panel3$statefip,FUN = mean)
panel3$female_s <- ave(panel3$female,panel3$statefip,FUN = mean)
panel3$yrseduc_s <- ave(panel3$yrseduc,panel3$statefip,FUN = mean)
panel3$lnwage_s <- ave(panel3$yrseduc,panel3$statefip,FUN= mean)

panel3$age_f <- panel3$age - panel3$age_s
panel3$age2_f <- panel3$age2 - panel3$age2_s
panel3$female_f <- panel3$female - panel3$female_s
panel3$yrseduc_f <- panel3$yrseduc - panel3$yrseduc_s
panel3$lnwage_f <- panel3$lnwage - panel3$lnwage_s


#####summaries of these variables###
#install.packages("psych")
library(psych)
panel3vars <- panel3[,c(10:19)]
describe(panel3vars)

par(mfrow=c(2,4))
plot(density(panel3$age_s))
plot(density(panel3$age2_s))
plot(density(panel3$female_s))
plot(density(panel3$yrseduc_s))
plot(density(panel3$age_f))
plot(density(panel3$age2_f))
plot(density(panel3$female_f))
plot(density(panel3$yrseduc_f))

##notes:
  #it looks like female_s has an outlier & female_f is skewed
  #

summary(panel3$age_s)
boxplot(panel3$age_s,panel3$female_s, panel3$yrseduc)

boxplot(panel3$age_f ~ panel3$statefip)
boxplot(panel3$age2_f ~ panel3$statefip)
boxplot(panel3$female_f ~ panel3$statefip)
boxplot(panel3$yrseduc_f ~ panel3$statefip)

