library(robumeta)
library(metafor)
library(clubSandwich)

# fit a simple meta-analysis model using robumeta

data(corrdat)
corrdat$studyid <- factor(corrdat$studyid)

robu_fit <- robu(effectsize ~ 1, data = corrdat, 
                 modelweights = "CORR", studynum = studyid,
                 var.eff.size = var)
robu_fit

# note value of tau^2
tau_sq_robu <- as.numeric(robu_fit$mod_info$tau.sq)

# fit the same specification using rma.uni
rma_uni_fit <- rma.uni(yi = effectsize, vi = var, data = corrdat)

# different estimate of tau^2
rma_uni_fit$tau2

# robust standard error from metafor
robust(rma_uni_fit, cluster = corrdat$studyid, adjust = TRUE)

# you can get the same robust SE from clubSandwich using vcov = "CR1"
coef_test(rma_uni_fit, cluster = corrdat$studyid, vcov = "CR1")

# Using the more accurate small-sample approximation vcov = "CR2"
# leads to slightly larger SE
coef_test(rma_uni_fit, cluster = corrdat$studyid, vcov = "CR2")

# calculate the inverse of the weight used by robumeta
corrdat$k <- with(corrdat, table(studyid)[studyid])
corrdat$Vbar <- with(corrdat, tapply(var, studyid, mean)[studyid])
corrdat$Vij <- with(corrdat, as.numeric(k * (Vbar + tau_sq_robu)))

# use the robumeta weights in rma.uni
rma_uni_robu <- rma.uni(yi = effectsize, vi = Vij, data = corrdat, method = "FE")

# average effect size estimates now agree
rma_uni_robu$beta
robu_fit$reg_table$b.r

# have to use "CR2" to get the standard errors to agree
coef_test(rma_uni_robu, cluster = corrdat$studyid, vcov = "CR2")
sqrt(vcovCR(rma_uni_robu, cluster = corrdat$studyid, type = "CR2"))
robu_fit$reg_table$SE

# the robumeta model is closer to the following multivariate model
V_mat <- with(corrdat, impute_covariance_matrix(vi = var, cluster = studyid, r = 0.7))
rma_mv_fit <- rma.mv(yi = effectsize, V = V_mat, data = corrdat,
                     random = ~ studyid)
coef_test(rma_mv_fit, cluster = corrdat$studyid, vcov = "CR2")
