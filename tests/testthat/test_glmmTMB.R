context("glmmTMB objects")

skip_if_not_installed("glmmTMB")
skip_if_not_installed("lme4")

data("sleepstudy", package = "lme4")

test_that("Gaussian random slopes agree with lmer", {
  tmb <- glmmTMB::glmmTMB(Reaction ~ Days + (Days | Subject), sleepstudy)
  lmer <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = FALSE)

  expect_equal(coef_CS(tmb), coef_CS(lmer), tolerance = 1e-4)
  expect_true(check_bread(tmb, sleepstudy$Subject, sleepstudy$Reaction))
  expect_true(check_CR(tmb, vcov = "CR2"))

  for (type in c("CR0", "CR2", "CR3")) {
    V <- vcovCR(tmb, type = type)
    expect_equal(as.matrix(V), as.matrix(vcovCR(lmer, type = type)), tolerance = 1e-4)
    expect_equal(as.matrix(vcovCR(tmb, cluster = sleepstudy$Subject, type = type)),
                 as.matrix(V))
    expect_equal(coef_test(tmb, vcov = V)$SE,
                 coef_test(lmer, vcov = type)$SE, tolerance = 1e-4)
  }
  expect_equal(Wald_test(tmb, constraints = constrain_zero(2), vcov = "CR2")$Fstat,
               Wald_test(lmer, constraints = constrain_zero(2), vcov = "CR2")$Fstat,
               tolerance = 1e-4)
  V <- vcovCR(tmb, type = "CR2")
  expect_equal(as.matrix(vcovCR(tmb, type = "CR2", target = attr(V, "target"))),
               as.matrix(V))
  expect_equal(unname(tcrossprod(vcovCR(tmb, type = "CR2", form = "estfun"))),
               unname(as.matrix(V)))
})

test_that("diagonal random slopes and REML agree with lmer", {
  tmb <- glmmTMB::glmmTMB(Reaction ~ Days + diag(1 + Days | Subject),
                          sleepstudy, REML = TRUE)
  lmer <- lme4::lmer(Reaction ~ Days + (Days || Subject), sleepstudy, REML = TRUE)
  expect_equal(as.matrix(vcovCR(tmb, type = "CR2")),
               as.matrix(vcovCR(lmer, type = "CR2")), tolerance = 1e-4)
})

test_that("dispersion, offsets and missing observations use the fitted data", {
  dat <- sleepstudy
  dat$off <- dat$Days / 2
  dat$Reaction[3] <- NA
  fit <- glmmTMB::glmmTMB(Reaction ~ Days + offset(off) + (1 | Subject),
                          dat, dispformula = ~ Days, na.action = na.exclude)
  frame <- model.frame(fit)
  rows <- split(seq_len(nrow(frame)), frame$Subject)
  residual_sd <- predict(fit, type = "disp")[-3]
  random_var <- as.numeric(glmmTMB::VarCorr(fit)$cond[[1]])
  expected <- lapply(rows, function(i) {
    matrix(random_var, length(i), length(i)) + diag(residual_sd[i]^2)
  })

  expect_equal(lapply(targetVariance(fit, frame$Subject), unname),
               lapply(expected, unname))
  expect_equal(residuals_CS(fit),
               as.numeric(frame$Reaction - frame$off -
                            model.matrix(fit) %*% glmmTMB::fixef(fit)$cond))
  expect_equal(as.matrix(vcovCR(fit, cluster = dat$Subject, type = "CR2")),
               as.matrix(vcovCR(fit, cluster = frame$Subject, type = "CR2")))
  expect_true(check_bread(fit, frame$Subject, frame$Reaction - frame$off))
  expect_true(check_CR(fit, vcov = "CR2"))
})

test_that("nested and crossed groups are checked against the clusters", {
  dat <- sleepstudy
  dat$school <- factor(rep(seq_len(6), each = 30))
  dat$crossed <- factor(rep(seq_len(10), 18))
  dat$all <- factor(rep(c("a", "b"), each = 90))
  nested <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | school) + (1 | Subject), dat)
  expect_equal(as.matrix(vcovCR(nested, type = "CR0")),
               as.matrix(vcovCR(nested, cluster = dat$school, type = "CR0")))
  expect_true(check_bread(nested, dat$school, dat$Reaction))
  expect_error(vcovCR(nested, cluster = dat$Subject, type = "CR0"), "nested")

  set.seed(467)
  dat$Reaction <- 100 + 2 * dat$Days + rnorm(18, sd = 15)[dat$Subject] +
    rnorm(10, sd = 12)[dat$crossed] + rnorm(nrow(dat), sd = 5)
  crossed <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject) + (1 | crossed), dat)
  expect_error(vcovCR(crossed, type = "CR0"), "outermost")
  expect_error(vcovCR(crossed, cluster = dat$all, type = "CR0"), "nested")
  dat$outer <- factor(rep(seq_len(2), each = 90))
  dat$crossed <- interaction(dat$outer, dat$crossed, drop = TRUE)
  crossed <- update(crossed, data = dat)
  V <- vcovCR(crossed, cluster = dat$outer, type = "CR0")
  expect_true(check_bread(crossed, dat$outer, dat$Reaction))
  expect_true(check_CR(crossed, vcov = vcovCR(crossed, cluster = dat$outer, type = "CR2")))
  expect_s3_class(V, "vcovCR")
})

test_that("models without random effects and rank-deficient fits work", {
  dat <- sleepstudy
  dat$duplicate <- dat$Days
  fit <- glmmTMB::glmmTMB(Reaction ~ Days + duplicate, dat,
                          control = glmmTMB::glmmTMBControl(rank_check = "adjust"))
  expect_error(vcovCR(fit, type = "CR0"), "clustering variable")
  expect_true(anyNA(coef_CS(fit)))
  lm_fit <- lm(Reaction ~ Days, dat)
  expect_equal(as.matrix(vcovCR(fit, cluster = dat$Subject, type = "CR2")),
               as.matrix(vcovCR(lm_fit, cluster = dat$Subject, type = "CR2")),
               tolerance = 1e-5)
})

test_that("unsupported fits fail explicitly", {
  dat <- sleepstudy
  dat$count <- as.integer(dat$Reaction / 100)
  poisson <- glmmTMB::glmmTMB(count ~ Days + (1 | Subject), dat, family = poisson)
  expect_error(vcovCR(poisson, type = "CR0"), "Gaussian")
  weighted <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject), dat,
                               weights = rep(1, nrow(dat)))
  expect_error(vcovCR(weighted, type = "CR0"), "weights")
  zi <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject), dat, ziformula = ~1)
  expect_error(vcovCR(zi, type = "CR0"), "Zero-inflated")
  mapped <- glmmTMB::glmmTMB(Reaction ~ Days + (1 | Subject), dat,
                             map = list(beta = factor(c(NA, 1))))
  expect_error(vcovCR(mapped, type = "CR0"), "mapped")
  ar <- glmmTMB::glmmTMB(Reaction ~ Days + ar1(factor(Days) + 0 | Subject), dat,
                         control = glmmTMB::glmmTMBControl(conv_check = "skip",
                                                           eigval_check = FALSE))
  expect_error(vcovCR(ar, type = "CR0"), "covariance structures")
})
