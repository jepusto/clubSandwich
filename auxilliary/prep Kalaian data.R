SATcoaching <- read.delim("auxilliary/Kalaian-Raudenbush-1996.txt")
SATcoaching$Outcome <- ifelse(SATcoaching$x==1, "SATM", "SATV")
SATcoaching$x <- NULL

save(SATcoaching, file = "data/SATcoaching.RData", compress = "xz")
head(SATcoaching)
dim(SATcoaching)
