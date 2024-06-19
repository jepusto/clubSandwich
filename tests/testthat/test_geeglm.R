context("geeglm objects")
set.seed(202201030)

skip_if_not_installed("geepack")

library(geepack)

J <- 20
tp <- 5
# Simulating a dataset
idvar <- rep(1:J, each=tp)
idvar2 <- factor(sample(LETTERS)[idvar])
timeorder <- rep(1:tp, J)
tvar      <- timeorder + rnorm(length(timeorder))
x1 <- rnorm(length(timeorder))
x2 <- rnorm(J)[idvar]
uuu   <- rep(rnorm(J), each=tp)
yvar  <- 1 + 0.4 * x1 + 0.2 * x2 + 2 * tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, x1, x2, tvar, yvar, idvar2)

simdatPerm <- simdat[sample(nrow(simdat)),]
simdatPerm <- simdatPerm[order(simdatPerm$idvar),]
wav <- simdatPerm$timeorder

# AR1 + wave
geeglm_AR1_wav <- geeglm(yvar ~ tvar + x1, id = idvar, 
                      data = simdatPerm, 
                      corstr = "ar1", waves = timeorder)
geeglm_AR1_wav

# AR1
geeglm_AR1 <- geeglm(yvar ~ tvar, id = idvar, data = simdat, family = "gaussian", corstr = "ar1")

# Independence
geeglm_ind <- geeglm(yvar ~ tvar + x1 + x2, id = idvar, data = simdat, corstr = "independence")

# Exchangeable
geeglm_exch <- geeglm(yvar ~ tvar + x2, id = idvar, data = simdat, corstr = "exchangeable")
geeglm_exch2 <- geeglm(yvar ~ tvar + x2, id = idvar2, data = simdat, corstr = "exchangeable")

# Unstructured
geeglm_unstr <- geeglm(yvar ~ tvar, id = idvar, data = simdat, corstr = "unstructured")
geeglm_unstr2 <- geeglm(yvar ~ tvar, id = idvar2, data = simdat, corstr = "unstructured")

# User-defined
zcor <- genZcor(clusz = table(simdat$idvar), 
                waves = simdat$timeorder, 
                corstrv = 4)
zcor_user <- zcor

geeglm_user <- geeglm(yvar ~ tvar + x1, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor_user, corstr = "userdefined")

# User-defined, Toeplitz
zcor_toep     <- matrix(NA, nrow(zcor), 4)
zcor_toep[,1] <- apply(zcor[,c(1, 5, 8,10)], 1, sum)
zcor_toep[,2] <- apply(zcor[,c(2, 6, 9)], 1, sum)
zcor_toep[,3] <- apply(zcor[,c(3, 7)], 1, sum)
zcor_toep[,4] <- zcor[,4]

geeglm_toep <- geeglm(yvar ~ tvar + x1, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor_toep, corstr = "userdefined")

# Fixed correlation
cor_fix <- matrix(c(1    , 0.5  , 0.25,  0.125, 0.125,
                    0.5  , 1    , 0.25,  0.125, 0.125,
                    0.25 , 0.25 , 1   ,  0.5  , 0.125,
                    0.125, 0.125, 0.5  , 1    , 0.125,
                    0.125, 0.125, 0.125, 0.125, 1     ), nrow=5, ncol=5)
zcor_fix <- fixed2Zcor(cor_fix, id=simdat$idvar, waves=simdat$timeorder)

geeglm_fix <- geeglm(yvar ~ tvar + x1 + x2, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor_fix, corstr = "fixed")

test_that("bread works", {
  
  expect_true(check_bread(geeglm_AR1_wav, cluster = simdatPerm$idvar, y = simdatPerm$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_AR1, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_ind, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_exch, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_exch2, cluster = simdat$idvar2, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_unstr, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_unstr2, cluster = simdat$idvar2, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_user, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_toep, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_fix, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  
})


test_that("vcovCR options work for CR2", {
  
  CR2_AR1_wav <- vcovCR(geeglm_AR1_wav, type = "CR2")
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR2"), CR2_AR1_wav)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR2", inverse_var = TRUE), CR2_AR1_wav)
  expect_false(identical(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = simdatPerm$idvar, inverse_var = FALSE), CR2_AR1_wav))
    
  target <- targetVariance(geeglm_AR1_wav, cluster = simdatPerm$idvar) 
  expect_equal(vcovCR(geeglm_AR1_wav, type = "CR2", target = target, inverse_var = TRUE), 
               CR2_AR1_wav, ignore_attr = TRUE)
  attr(CR2_AR1_wav, "inverse_var") <- FALSE
  expect_equal(vcovCR(geeglm_AR1_wav, type = "CR2", target = target, inverse_var = FALSE), 
               CR2_AR1_wav, ignore_attr = TRUE)
    
})
  
  
test_that("vcovCR options work for CR4", {
  CR4_AR1_wav <- vcovCR(geeglm_AR1_wav, type = "CR4")
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4"), CR4_AR1_wav)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = TRUE), CR4_AR1_wav)
  expect_false(identical(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = FALSE), CR4_AR1_wav))
  
  target <- targetVariance(geeglm_AR1_wav, cluster = simdatPerm$idvar)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", 
                      target = target, inverse_var = TRUE), CR4_AR1_wav, ignore_attr = TRUE)
  attr(CR4_AR1_wav, "inverse_var") <- FALSE
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", 
                      target = target, inverse_var = FALSE), CR4_AR1_wav, ignore_attr = TRUE)
})


test_that("CR2 and CR4 are target-unbiased", {
  
  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_AR1, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_ind, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_exch, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_exch2, cluster = simdat$idvar2, vcov = "CR2"))
  expect_true(check_CR(geeglm_unstr, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_unstr2, cluster = simdat$idvar2, vcov = "CR2"))
  expect_true(check_CR(geeglm_user, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_toep, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_fix, cluster = simdat$idvar, vcov = "CR2"))

  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_AR1, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_ind, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_exch, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_exch2, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_unstr, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_unstr2, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_user, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_toep, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_fix, cluster = simdat$idvar, vcov = "CR4"))
  
})



CR_types <- paste0("CR",0:4)

test_that("Order doesn't matter.", {
  
  skip_on_cran()
  
  check_sort_order(geeglm_ind, dat = simdat, CR_types = CR_types)
  check_sort_order(geeglm_ind, dat = simdat, arrange = "idvar2", CR_types = CR_types)
  check_sort_order(geeglm_exch, dat = simdat, arrange = "idvar", CR_types = CR_types)
  check_sort_order(geeglm_exch2, dat = simdat, arrange = "idvar", CR_types = CR_types)
  check_sort_order(geeglm_exch2, dat = simdat, arrange = "idvar2", CR_types = CR_types)
  expect_error(check_sort_order(geeglm_unstr, dat = simdat, arrange = "idvar", CR_types = CR_types))

})

  
test_that("clubSandwich works with dropped observations", {
  
  dat_miss <- simdat
  dat_miss$yvar[sample.int(nrow(simdat), size = round(nrow(simdat) / 10))] <- NA
  dat_complete <- subset(dat_miss, !is.na(yvar))
  
  mod_dropped <- geeglm(yvar ~ tvar, id = idvar, 
                       data = dat_miss, corstr = "independence")
  mod_complete <- geeglm(yvar ~ tvar, id = idvar, 
                                 data = dat_complete, corstr = "independence")
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(mod_dropped, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(mod_complete, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(mod_dropped, cluster = dat_miss$idvar, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(mod_complete, cluster = dat_complete$idvar, vcov = x, test = "All", p_values = FALSE))
  compare_ttests(test_drop, test_complete)
  
})


test_that("vcovCR works for clustering variables higher than id variable.", {
  
  # create higher level
  pair_id <- rep(1:nlevels(as.factor(simdat$idvar)), each = 3, length.out = nlevels(as.factor(simdat$idvar)))[as.factor(simdat$idvar)] # factor cluster
  pair_id <- factor(pair_id)
  
  # cluster at higher level
  V_AR1 <- vcovCR(geeglm_AR1, type = "CR2", cluster =  pair_id)
  V_AR1_wav <- vcovCR(geeglm_AR1_wav, type = "CR2", cluster =  pair_id)
  V_ind <- vcovCR(geeglm_ind, type = "CR2", cluster =  pair_id)
  V_exch <- vcovCR(geeglm_exch, type = "CR2", cluster =  pair_id)
  V_exch2 <- vcovCR(geeglm_exch2, type = "CR2", cluster =  pair_id)
  V_unstr <- vcovCR(geeglm_unstr, type = "CR2", cluster =  pair_id)
  V_unstr2 <- vcovCR(geeglm_unstr2, type = "CR2", cluster =  pair_id)
  V_user <- vcovCR(geeglm_user, type = "CR2", cluster =  pair_id)
  V_toep <- vcovCR(geeglm_toep, type = "CR2", cluster =  pair_id)
  V_fix <- vcovCR(geeglm_fix, type = "CR2", cluster =  pair_id)
  expect_is(V_AR1, "vcovCR")
  expect_is(V_AR1_wav, "vcovCR")
  expect_is(V_ind, "vcovCR")
  expect_is(V_exch, "vcovCR")
  expect_is(V_exch2, "vcovCR")
  expect_is(V_unstr, "vcovCR")
  expect_is(V_unstr2, "vcovCR")
  expect_is(V_user, "vcovCR")
  expect_is(V_toep, "vcovCR")
  expect_is(V_fix, "vcovCR")
  
  # check that result does not depend on sort-order
  scramble_id <- factor(simdat$idvar, levels = sample(1:J))
  dat_higher <- cbind(simdat, pair_id = factor(pair_id), scramble_id)
  dat_higher <- dat_higher[order(dat_higher$scramble_id),]
  
  dat_scramble <- dat_higher[sample(nrow(dat_higher)),]
  dat_scramble <- dat_scramble[order(dat_scramble$scramble_id),]

  V_AR1_scramble <- vcovCR(update(geeglm_AR1, data = dat_higher), 
                               type = "CR2", cluster = dat_higher$pair_id)
  expect_equal(V_AR1, V_AR1_scramble, tol = 10^-6, check.attributes = FALSE)
  
  V_AR1_wav_scramble <- vcovCR(update(geeglm_AR1_wav, data = dat_higher), 
                       type = "CR2", cluster = dat_higher$pair_id)
  expect_equal(V_AR1_wav, V_AR1_wav_scramble, tol = 10^-6, check.attributes = FALSE)
  
  V_ind_scramble <- vcovCR(update(geeglm_ind, data = dat_scramble), 
                               type = "CR2", cluster = dat_scramble$pair_id)
  expect_equal(V_ind, V_ind_scramble, tol = 10^-6, check.attributes = FALSE)
  
  V_exch2_scramble <- vcovCR(update(geeglm_exch2, data = dat_scramble), 
                            type = "CR2", cluster = dat_scramble$pair_id)
  expect_equal(V_exch2, V_exch2_scramble, tol = 10^-6, check.attributes = FALSE)
  
  V_unstr_scramble <- vcovCR(update(geeglm_unstr, data = dat_higher), 
                           type = "CR2", cluster = dat_higher$pair_id)
  expect_equal(V_unstr, V_unstr_scramble, tol = 10^-6, check.attributes = FALSE)
  
})

check_geeglm <- function(obj) {
  
  cr0_pack <- vcov(obj)
  cr0_club <- as.matrix(vcovCR(obj, type = "CR0"))
  
  expect_equal(cr0_pack, cr0_club)
  
  cr3_pack <- vcov(update(obj, std.err = "jack"))
  cr3_club <- as.matrix(vcovCR(obj, type = "CR3"))
  
  J <- length(unique(obj$id))
  p <- nrow(obj$geese$infls)
  f <- (J - p) / J

  expect_equal(cr3_pack, f * cr3_club)
}

test_that("vcovCR agrees with geeglm for CR0 and CR3.", {
  
  check_geeglm(geeglm_AR1_wav)
  check_geeglm(geeglm_AR1)
  check_geeglm(geeglm_ind)
  check_geeglm(geeglm_exch)
  check_geeglm(geeglm_exch2)
  check_geeglm(geeglm_unstr)
  check_geeglm(geeglm_unstr2)
  check_geeglm(geeglm_user)
  check_geeglm(geeglm_toep)
  check_geeglm(geeglm_fix)

})
