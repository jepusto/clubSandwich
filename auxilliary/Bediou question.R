load("auxilliary/es_data.RData")

## Multilevel model
vcov_mat <- impute_covariance_matrix(es_data$g_var, 
                                     cluster = es_data$Paper,                            
                                     r = 0.8) 
mlma_model<- metafor::rma.mv(g ~ 1, # Intercept only MLMA model
                             V = vcov_mat, 
                             random = ~ 1 | Paper / ES_ID,
                             sparse = TRUE,
                             data = es_data, 
                             test="t")
## Add RVE using clubSandwich
library(clubSandwich)

CHE_model <- clubSandwich::coef_test(mlma_model, vcov = "CR2")  