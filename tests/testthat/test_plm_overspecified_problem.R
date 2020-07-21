context("plm objects - balanced panel with cluster-level interactions")

library(plm, quietly=TRUE)

set.seed(20200721)
G <- 93
N <- 4 * G
Ts <- 4

beta1 <- rgamma(Ts, shape = 0.3, rate = 0.1)
beta2 <- rgamma(Ts, shape = 0.2, rate = 0.1)

grp <- factor(rep(1:G, each = 4 * Ts))
ID <- factor(rep(1:N, each = Ts))
trt <- factor(rep(LETTERS[1:4], times = N))
X1 <- rep(rnorm(N), each = Ts)
X2 <- rep(rnorm(N), each = Ts)
Y <- beta1[trt] * X1 + beta2[trt] * X2 + rnorm(N)[ID] + rnorm(N * Ts, sd = 0.5) 

dat <- data.frame(grp, ID, trt, X1, X2, Y)
dat <- pdata.frame(dat, index = c("ID","trt"))

obj <- plm(Y ~ trt + trt * X1 + trt * X2, 
           data=dat, 
           model="within")

# cluster <- dat$grp
# type <- "CR2"
# target <- NULL
# inverse_var <- TRUE
# form <- "sandwich"
# ignore_FE <- FALSE
# 
# 
# colnames(model_matrix.plm(obj))
# model.matrix(obj) %>% str()
# model.matrix(obj, cstcovar.rm = "all") %>% str()

CR_types <- paste0("CR",0:2)

test_that("vcovCR works with cluster-level interactions.", {
  meat_list <- lapply(CR_types, function(x) vcovCR(obj = obj, cluster=dat$grp, type = x, form = "meat"))
  bread_dim <- dim(bread(obj))
  lapply(meat_list, function(x) expect_identical(dim(x), bread_dim))
  
  V_CR_list <- lapply(CR_types, function(x) vcovCR(obj = obj, cluster=dat$grp, type = x))
  lapply(V_CR_list, expect_s3_class, class = "vcovCR")
  
})

test_that("CR0 agrees with built-in CRVE for plm", {
  V_plm <- vcovHC(obj, method = "arellano", type = "HC0")
  V_CR0 <- vcovCR(obj, type = "CR0")
  expect_equal(V_plm, as.matrix(V_CR0), check.attributes = FALSE)
})
