library(dplyr)
library(haven)
library(plm)
library(clubSandwich)

Shores <-
  read_dta("auxilliary/Shores_test.dta") %>%
  rename(W = `_W_Weight`)

Shores_fit <- plm(delta ~ treat + factor(year), data = Shores,
                  weights = W,
                  index = c("fips","year"),
                  effect = "individual")
summary(Shores_fit)
coef_test(Shores_fit, vcov = "CR2")
