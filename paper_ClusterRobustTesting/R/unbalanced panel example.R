library(foreign)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plm)

deaths <- read.dta("paper_ClusterRobustTesting/data/deaths.dta") 

filter(deaths, agegr=="18-20 yrs" & year <= 1983 & dtype == "MVA" & !is.na(beertaxa)) %>%
  select(-agegr, -dtype) %>%
  group_by(state) %>%
  mutate(legal_s = legal - mean(legal),
         beertaxa_s = beertaxa - mean(beertaxa)) %>%
  as.data.frame() ->
  death_dat

death_dat$mrate_miss <- death_dat$mrate
death_dat$mrate_miss[sample(nrow(death_dat), size = 100)] <- NA

var_names <- c("legal","beertaxa")
lm_fit <- lm(mrate_miss ~ 0 + legal + beertaxa + factor(state) + factor(year), data = death_dat)
CR2_iv_lm <- vcovCR(lm_fit, cluster = death_dat$state, type = "CR2")
CR2_iv_lm[var_names, var_names]
vcovCR(lm_fit, cluster = death_dat$state, type = "CR2", inverse_var = FALSE)[var_names, var_names]
plm_fit <- plm(mrate_miss ~ legal + beertaxa, data = death_dat, 
               index = c("state", "year"), effect = "twoways", model = "within")
coef(lm_fit)[var_names]
coef(plm_fit)
vcovCR(plm_fit, cluster = "individual", type = "CR2")
vcovCR(plm_fit, cluster = death_dat$state, type = "CR2", inverse_var = FALSE)


#------------------------
# Compute A matrices
#------------------------

cluster <- with(death_dat, as.factor(state[!is.na(mrate_miss)]))

X <- model.matrix(lm_fit)
X_names <- colnames(X)
R <- X[,var_names]
S <- X[,grepl("year",X_names)]
T <- X[,grepl("state",X_names)]

Sp <- residuals(lm.fit(T, S))
Rp <- residuals(lm.fit(Sp, residuals(lm.fit(T, R))))
Up <- residuals(lm.fit(T, cbind(R, S)))

R_list <- matrix_list(Rp, cluster, "row")
U_list <- matrix_list(Up, cluster, "row")
T_list <- matrix_list(T, cluster, "row")
M_R <- chol2inv(chol(crossprod(Rp)))
M_U <- chol2inv(chol(crossprod(Up)))
M_T <- chol2inv(chol(crossprod(T)))
IH_list <- IH_jj_list(M = M_U, X_list = U_list, XW_list = lapply(U_list, t))
A_mats <- lapply(IH_list, Sym_power, p = -1/2)
terms <- Map(function(t_mat, a, r) t(t_mat) %*% a %*% r, t_mat = T_list, a = A_mats, r = R_list)
meat <- Reduce("+", lapply(terms, function(x) t(x) %*% M_T %*% x))
M_R %*% meat %*% M_R
