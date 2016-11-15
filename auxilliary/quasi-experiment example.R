library(plm)
devtools::load_all()
rm(list=ls())
load("auxilliary/quasi-experiment.Rdat")

d_p <- plm.data(d, indexes=c("sgrp"))

p <- plm(y ~ treat + frl, data = d_p, effect="individual", model="within", index="sgrp")

system.time(coef_test(p, vcov = "CR1", cluster = d$sgrp, test = "naive-t"))
system.time(coef_test(p, vcov = "CR2", cluster = d$sgrp, test = "naive-t"))
system.time(coef_test(p, vcov = "CR1", cluster = d$sgrp, test = "Satterthwaite"))
system.time(coef_test(p, vcov = "CR2", cluster = d$sgrp, test = "Satterthwaite"))
# system.time(coef_test_old(p, vcov = "CR1", cluster = d$sgrp, test = "Satterthwaite"))
# system.time(coef_test_old(p, vcov = "CR2", cluster = d$sgrp, test = "Satterthwaite"))
