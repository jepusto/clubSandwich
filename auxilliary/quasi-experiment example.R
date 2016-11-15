library(plm)
library(clubSandwich)
# devtools::load_all()
rm(list=ls())
load("auxilliary/quasi-experiment.Rdata")

d_p <- plm.data(d, indexes=c("sgrp"))

p <- plm(y ~ treat + frl, data = d_p, effect="individual", model="within", index="sgrp")

system.time(CR1_t <- coef_test(p, vcov = "CR1", cluster = d$sgrp, test = "naive-t"))
system.time(CR2_t <- coef_test(p, vcov = "CR2", cluster = d$sgrp, test = "naive-t"))
system.time(CR1_Satt <- coef_test(p, vcov = "CR1", cluster = d$sgrp, test = "Satterthwaite"))
system.time(CR2_Satt <- coef_test(p, vcov = "CR2", cluster = d$sgrp, test = "Satterthwaite"))
# system.time(CR1_Satt_old <- coef_test_old(p, vcov = "CR1", cluster = d$sgrp, test = "Satterthwaite"))
# system.time(CR2_Satt_old <- coef_test_old(p, vcov = "CR2", cluster = d$sgrp, test = "Satterthwaite"))
