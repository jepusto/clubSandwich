context("gls objects")

library(nlme)

data(Ovary, package = "nlme")

fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
