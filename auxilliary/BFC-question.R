library(metafor)
library(multcomp)
library(clubSandwich)

harm <- read.table(file="auxilliary/harm.txt")

meta_res <-rma.mv(y=r, V=var, mods= ~ sampletype1 - 1, random = ~ 1 | study/outcome, data=harm)

rob_er_meta <- conf_int(meta_res, vcov="CR2", level = 0.95, test = "Satterthwaite")

summary(glht(meta_res, linfct=cbind(contrMat(rep(1,4), type="Tukey"))))

constraint_mat <- cbind(rep(-1, 3), diag(1, 3))
Wald_test(meta_res, constraints = constraint_mat, vcov = "CR2")
