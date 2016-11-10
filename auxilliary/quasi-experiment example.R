library(plm)
devtools::load_all()

rm(list=ls())
set.seed(20161109)

Nt <- 500
Cper <- 50
N <- (1 + Cper) * Nt


dat <- data.frame(
  sid = 1:51,
  trt = rep(c(1, rep(0,Cper)), Nt),
  matched = rep(1:Nt, each = Cper + 1),
  school = sample(5:1000, size = N, replace = TRUE)
)

dat <- within(dat, {
  school[trt==1] <- sample(1:4, size = Nt, replace = TRUE)
  matched <- factor(matched)
  school <- factor(school)
  y <- trt + rnorm(nlevels(school))[school] + rnorm(N)
})

obj <- plm(y ~ trt, data = dat, effect = "individual", index = c("matched","sid"))
summary(obj)
vcovHC(obj, method = "arellano", type = "HC0")
vcovCR(obj, cluster = dat$match, type = "CR0")

coef_test(obj, vcov = "CR0", cluster = dat$school, test = "z")
coef_test(obj, vcov = "CR2", cluster = dat$school, test = "z")
system.time(cr0_school <- coef_test(obj, vcov = "CR0", cluster = dat$school, test = "Satterthwaite"))

vcr <- vcovCR(obj, cluster = dat$school, type = "CR2")

cluster <- dat$school
type <- "CR2"
target <- NULL
inverse_var <- is.null(target)
form <- "sandwich"
obj$na.action <- attr(obj$model, "na.action")
vcov <- vcov_CR(obj, cluster = cluster, 
                type = type, 
                target = target, 
                inverse_var = inverse_var, 
                form = form)
test <- "Satterthwaite"
