
impute_covariance_matrix <- function(vi, cluster, r, return_list = identical(as.factor(cluster), sort(as.factor(cluster)))) {
  
  vi_list <- split(vi, cluster)
  r_list <- rep_len(r, length(vi_list))
  vcov_list <- Map(function(V, rho) (rho + diag(1 - rho, nrow = length(V))) * tcrossprod(sqrt(V)), V = vi_list, rho = r_list)
  
  if (return_list) {
    return(vcov_list)
  } else {
    vcov_mat <- metafor::bldiag(vcov_list)
    cluster_index <- order(order(cluster))
    return(vcov_mat[cluster_index, cluster_index])
  }
}

library(dplyr)
data(SATcoaching, package = "clubSandwich")

mean_hrs_ln <- 
  SATcoaching %>% 
  group_by(study) %>%
  summarise(hrs_ln = mean(log(hrs))) %>%
  summarise(hrs_ln = mean(hrs_ln, na.rm = TRUE))

SATcoaching <- 
  SATcoaching %>%
  # clean variables
  mutate(
    study = as.factor(study),
    hrs_ln = log(hrs) - mean_hrs_ln$hrs_ln
  ) %>%
  # sort by study ID
  arrange(study, test)

V_list <- with(SATcoaching, impute_covariance_matrix(vi = V, cluster = study, r = 0.66))
MVFE_hrs <- rma.mv(d ~ 0 + test + test:hrs_ln, V = V_list, 
                   data = SATcoaching)
coef_test(MVFE_hrs, vcov = "CR2", cluster = SATcoaching$study)

lm_fit <- lm(d ~ 0 + test + test:hrs_ln, data = SATcoaching)
coef_test(lm_fit, vcov = "CR2", cluster = SATcoaching$study)
rma_fit <- rma.uni(d ~ 0 + test + test:hrs_ln, data = SATcoaching, vi = V)
coef_test(rma_fit, vcov = "CR2", cluster = SATcoaching$study)


MVRE_hrs <- rma.mv(d ~ 0 + test + test:hrs_ln, V = V_list, 
                   data = SATcoaching,
                   random = ~ test | study, struct = "UN")

V_mat <- with(SATcoaching, impute_covariance_matrix(vi = V, cluster = study, r = 0.66, return_list = FALSE))
has_hrs <- !is.na(SATcoaching$hrs_ln)
V_mat <- V_mat[has_hrs, has_hrs]
MVRE_hrs_filtered <- rma.mv(d ~ 0 + test + test:hrs_ln, V = V_mat, 
                            data = filter(SATcoaching, !is.na(hrs_ln)),
                            random = ~ test | study, struct = "UN")
coef_test(MVRE_hrs, vcov = "CR2")
coef_test(MVRE_hrs_filtered, vcov = "CR2")
