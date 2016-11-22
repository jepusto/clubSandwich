library(plm)
devtools::load_all()
rm(list=ls())

trt_size <- c(547, 709, 385, 129, 302, 537, 747, 218, 73, 373, 440, 78, 272, 98, 67, 211, 19, 39, 89, 24, 217, 79, 83, 157, 118)
FRL <- c(82, 51, 41, 83, 82, 60, 46, 49, 57, 90, 92, 99, 81, 90, 87, 100, 58, 58, 34, 61, 77, 73, 38, 91, 90)

(Ntrt_schools <- length(trt_size))
(Ntrt <- sum(trt_size))
(Nschools <- 5889 + Ntrt_schools)
Cper <- 50
(N <- (1 + Cper) * Ntrt)


dat <- data.frame(
  sid = 1:51,
  trt = rep(c(1, rep(0,Cper)), Ntrt),
  matched = rep(1:Ntrt, each = Cper + 1),
  school = sample(5:Nschools, size = N, replace = TRUE)
)

dat <- within(dat, {
  school[trt==1] <- rep(1:Ntrt_schools, trt_size)
  matched <- factor(matched)
  school <- factor(school)
  x <- rnorm(N)
  sex <- sample(0:1, size = N, replace = TRUE)
  y <- trt + rnorm(Nschools)[school] + rnorm(N)
})

table(table(dat$school))

obj <- plm(y ~ trt, data = dat, effect = "individual", model = "pooling", index = c("matched","sid"))
obj <- lm(y ~ trt, data = dat)
summary(obj)
CR2_t <- coef_test(obj, vcov = "CR2", cluster = dat$school, test = "Satterthwaite", ignore_FE = TRUE)
CR2_t
CR2_t$df
CR0_t <- coef_test(obj, vcov = "CR0", cluster = dat$school, test = "Satterthwaite")

NT <- sum(trt_size)
fi <- trt_size / NT
1 / sum(fi^2)
1 / (sum(fi^2) + sum(tcrossprod(fi^2 / (1 - fi))))


ctl_size <- with(dat, as.numeric(table(school[trt==0])))
NC <- sum(ctl_size)
gi <- ctl_size / NC
1 / sum(gi^2)
(1 / NT + 1 / NC)^2 / 
  (sum(fi^2) / NT^2 + sum(tcrossprod(fi^2 / (1 - fi))) / NT^2 + 
     sum(gi^2) / NC^2 + sum(tcrossprod( gi^2 / (1 - gi))) / NC^2)

#------------------------
# smaller example
#------------------------
rm(list=ls())
set.seed(20161109)

Ntrt_schools <- 4
Nschools <- 100
Ntrt <- 30
Cper <- 50
N <- (1 + Cper) * Ntrt


dat <- data.frame(
  sid = 1:51,
  trt = rep(c(1, rep(0,Cper)), Ntrt),
  matched = rep(1:Ntrt, each = Cper + 1),
  school = sample(5:Nschools, size = N, replace = TRUE)
)

dat <- within(dat, {
  school[trt==1] <- sample(1:Ntrt_schools, size = Ntrt, replace = TRUE)
  matched <- factor(matched)
  school <- factor(school)
  x <- rnorm(N)
  sex <- sample(0:1, size = N, replace = TRUE)
  y <- trt + rnorm(Nschools)[school] + rnorm(N)
})
table(table(dat$school))

obj <- plm(y ~ trt + x, data = dat, effect = "individual", index = c("matched","sid"))

summary(obj)
system.time(CR2_t <- coef_test(obj, vcov = "CR2", cluster = dat$school, test = "naive-t"))
system.time(CR2a_Satt <- coef_test(obj, vcov = "CR2", cluster = dat$school, test = "Satterthwaite"))
system.time(CR2b_Satt <- coef_test(obj, vcov = "CR2", cluster = dat$school, test = "Satterthwaite", ignore_FE = TRUE))
CR2_t
CR2a_Satt
CR2b_Satt

cluster <- dat$school
type <- "CR2"
target <- NULL
inverse_var <- is.null(target)
form <- "sandwich"
test <- "Satterthwaite"

obj$na.action <- attr(obj$model, "na.action")
vcov <- vcov_CR(obj, cluster = cluster, 
                type = type, 
                target = target, 
                inverse_var = inverse_var, 
                form = form)

beta <- coef_CS(obj)
beta_NA <- is.na(beta)
SE <- sqrt(diag(vcov))

system.time(S_array <- get_S_array(obj, vcov))
system.time(P_array <- get_P_array(GH <- get_GH(obj, vcov)))

Satterthwaite_old(beta, SE, S_array)
Satterthwaite(beta, SE, P_array)

saddlepoint_old(t_stats = beta[!beta_NA] / SE, S_array = S_array)
saddlepoint(t_stats = beta[!beta_NA] / SE, P_array = P_array)


constraints <- 1:2
test <- c("chi-sq","Naive-F","HTA","HTB","HTZ","EDF","EDT")
W_new <- Wald_test(obj, constraints = constraints, vcov = vcov, test = test)
W_old <- Wald_test_old(obj, constraints = constraints, vcov = vcov, test = test)

C_mat <- get_constraint_mat(obj, constraints)
GH <- get_GH(obj, vcov)
P_array <- get_P_array(GH, all_terms = TRUE)

q <- nrow(C_mat)
N <- dim(S_array)[2]
J <- dim(S_array)[3]
S_array <- array(apply(S_array, 3, function(s) C_mat %*% s), dim = c(q, N, J))
Omega <- apply(array(apply(S_array, 3, tcrossprod), dim = c(q,q,J)), 1:2, sum)
Omega_nsqrt <- matrix_power(Omega, -1/2)

Cov_arr_old <- covariance_array_old(S_array, Omega_nsqrt, q = q, J = J)
Cov_arr <- covariance_array(P_array, Omega_nsqrt)
identical(Cov_arr, Cov_arr_old)
all.equal(Cov_arr, Cov_arr_old)
Cov_arr / Cov_arr_old

tot_var_old <- total_variance_mat_old(S_array, Omega_nsqrt, q = q, J = J)
tot_var <- total_variance_mat(P_array, Omega_nsqrt)
identical(tot_var, tot_var_old)
all.equal(tot_var, tot_var_old)
tot_var / tot_var_old

Wald_testing(C_mat, beta, vcov, test, GH)
Wald_testing_old(C_mat, beta, vcov, test, S_array)


attr(vcov, "inverse_var") <- FALSE
system.time(S_array <- get_S_array(obj, vcov))
system.time(GH <- get_GH(obj, vcov))
system.time(P_array <- get_P_array(GH))
Satterthwaite_old(beta, SE, S_array)
Satterthwaite(beta, SE, P_array)
saddlepoint_old(t_stats = beta[!beta_NA] / SE, S_array = S_array)
saddlepoint(t_stats = beta[!beta_NA] / SE, P_array = P_array)

constraints <- 1:2
test <- c("chi-sq","Naive-F","HTA","HTB","HTZ","EDF","EDT")
W_new <- Wald_test(obj, constraints = constraints, vcov = vcov, test = test)
W_old <- Wald_test_old(obj, constraints = constraints, vcov = vcov, test = test)
W_new / W_old
