library(plm)
# library(clubSandwich)
devtools::load_all()

rm(list=ls())
load("auxilliary/quasi-experiment.Rdata")

d_p <- plm.data(d, indexes=c("sgrp"))

p <- plm(y ~ treat + frl, data = d_p, effect="individual", model="within", index="sgrp")

system.time(CR1_t <- coef_test(p, vcov = "CR1", cluster = d_p$cid, test = "naive-t"))
system.time(CR1a_Satt <- coef_test(p, vcov = "CR1", cluster = d_p$cid, test = "Satterthwaite"))
system.time(CR1b_Satt <- coef_test(p, vcov = "CR1", cluster = d_p$cid, test = "Satterthwaite", ignore_FE = TRUE))
# system.time(CR1_Satt_old <- coef_test_old(p, vcov = "CR1", cluster = d$cid, test = "Satterthwaite"))
CR1_t
CR1a_Satt
CR1b_Satt

system.time(CR2_t <- coef_test(p, vcov = "CR2", cluster = d_p$cid, test = "naive-t"))
system.time(CR2a_Satt <- coef_test(p, vcov = "CR2", cluster = d_p$cid, test = "Satterthwaite"))
system.time(CR2b_Satt <- coef_test(p, vcov = "CR2", cluster = d_p$cid, test = "Satterthwaite", ignore_FE = TRUE))
# system.time(CR2_Satt_old <- coef_test_old(p, vcov = "CR2", cluster = d$cid, test = "Satterthwaite"))
CR2_t
CR2a_Satt
CR2b_Satt

library(dplyr)
d %>% 
  group_by(cid) %>% 
  summarise(n = n(), trtd = mean(treat)) %>% 
  group_by(trtd) %>%
  summarise(schools = n(), students = sum(n))
