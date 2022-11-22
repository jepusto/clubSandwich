skip_if_not_installed("fixest")

library(fixest, quietly=TRUE)

set.seed(20221122)
data("airquality")
airquality <- airquality[complete.cases(airquality),]
airquality$wt <- 1L + rpois(nrow(airquality), lambda = 3)

fe_intercept <- feols(Ozone ~ Solar.R + Wind | Month, 
                      data = airquality)
fe_slope <- feols(Ozone ~ Solar.R + Wind | Month[Temp], 
                  data = airquality)
few_intercept <- feols(Ozone ~ Solar.R + Wind | Month, 
                       weights = ~ wt,
                       data = airquality)
few_slope <- feols(Ozone ~ Solar.R + Wind | Month[Temp], 
                   weights = ~ wt,
                   data = airquality)

model.matrix(fe_intercept)

data("base_did")
est_did <- feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
summary(est_did)
