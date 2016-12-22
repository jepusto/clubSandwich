library(dplyr)
library(tidyr)
library(readr)

SATcoaching <- 
  read_csv("auxilliary/Kalaian-Raudenbush-1996.csv") %>%
  gather("test","d", SATV, SATM) %>%
  mutate(test = ifelse(test == "SATV", "Verbal", "Math")) %>%
  filter(!is.na(d)) %>%
  arrange(study_type, study, test) %>%
  mutate(nT = ifelse(study=="Whitla" & test=="Math", 50, nT),
         nC = ifelse(study=="Whitla" & test=="Math", 50, nC),
         V = round(1 / nT + 1 / nC + d^2 / (2 * (nC + nT)), 4)) %>%
  select(study, year, test, d, V, nT, nC, study_type, hrs, ETS, homework) %>%
  as.data.frame()

save(SATcoaching, file = "data/SATcoaching.RData", compress = "xz")
head(SATcoaching)
dim(SATcoaching)
