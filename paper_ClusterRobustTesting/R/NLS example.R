setwd("paper_ClusterRobustTesting")
library(tidyr)
library(dplyr)
library(foreign)
library(plm)
nlswork <- read.dta("data/nlswork.dta", convert.factors=FALSE)
  