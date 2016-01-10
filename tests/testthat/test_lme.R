context("lme objects")
library(nlme)
data("Orthodont", package = "nlme")
fm1 <- lme(distance ~ age, random = ~ age, data = Orthodont)
fm2 <- lme(distance ~ age + Sex, random = ~ 1, data = Orthodont)
summary(fm1)
summary(fm2)

# test cluster specification
# test target matrix specification
# test inverse_var detection