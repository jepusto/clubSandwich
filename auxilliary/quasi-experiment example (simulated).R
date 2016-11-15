library(plm)
devtools::load_all()

rm(list=ls())
set.seed(20161109)

Nschools <- 100
Nt <- 30
Cper <- 50
N <- (1 + Cper) * Nt


dat <- data.frame(
  sid = 1:51,
  trt = rep(c(1, rep(0,Cper)), Nt),
  matched = rep(1:Nt, each = Cper + 1),
  school = sample(5:Nschools, size = N, replace = TRUE)
)

dat <- within(dat, {
  school[trt==1] <- sample(1:4, size = Nt, replace = TRUE)
  matched <- factor(matched)
  school <- factor(school)
  x <- rnorm(N)
  sex <- sample(0:1, size = N, replace = TRUE)
  y <- trt + rnorm(Nschools)[school] + rnorm(N)
})
table(table(dat$school))

obj <- plm(y ~ trt + x, data = dat, effect = "individual", index = c("matched","sid"))

summary(obj)
system.time(coef_test(obj, vcov = "CR1", cluster = dat$school, test = "naive-t"))
system.time(coef_test(obj, vcov = "CR2", cluster = dat$school, test = "naive-t"))
# system.time(coef_test_old(obj, vcov = "CR1", cluster = dat$school, test = "Satterthwaite"))
# system.time(coef_test_old(obj, vcov = "CR2", cluster = dat$school, test = "Satterthwaite"))
system.time(coef_test(obj, vcov = "CR1", cluster = dat$school, test = "Satterthwaite"))
system.time(coef_test(obj, vcov = "CR2", cluster = dat$school, test = "Satterthwaite"))


system.time(V_cr2 <- vcovCR(obj, cluster = dat$school, type = "CR2"))
system.time(coef_test(obj, vcov = V_cr2, test = "Satterthwaite", verbose = TRUE))


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
C_mat <- get_constraint_mat(obj, constraints)
test <- c("chi-sq","Naive-F","HTA","HTB","HTZ","EDF","EDT")
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

tot_var_old <- total_variance_mat_old(S_array, Omega_nsqrt, q = q, J = J)
tot_var <- total_variance_mat(P_array, Omega_nsqrt)
identical(tot_var, tot_var_old)
all.equal(tot_var, tot_var_old)
tot_var / tot_var_old

Wald_testing(C_mat, beta, vcov, test, P_array)


attr(vcov, "inverse_var") <- FALSE
system.time(S_array <- get_S_array(obj, vcov))
system.time(P_array <- get_P_array(obj, vcov, verbose = TRUE))
system.time(GH <- get_GH(obj, vcov))
Satterthwaite_old(beta, SE, P_array)
Satterthwaite(beta, SE, GH)
saddlepoint_old(t_stats = beta[!beta_NA] / SE, P_array = P_array)
saddlepoint(t_stats = beta[!beta_NA] / SE, GH = GH)


