context("geeglm objects")
set.seed(202201030)

library(geepack)

J <- 30
tp <- 5
# Simulating a dataset
timeorder <- rep(1:tp, J)
tvar      <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:J, each=tp)
uuu   <- rep(rnorm(J), each=tp)
yvar  <- 1 + 2 * tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)

simdatPerm <- simdat[sample(nrow(simdat)),]
simdatPerm <- simdatPerm[order(simdatPerm$idvar),]
wav <- simdatPerm$timeorder

# AR1 + wave
geeglm_AR1_wav <- geeglm(yvar ~ tvar, id = idvar, 
                      data = simdatPerm, 
                      corstr = "ar1", waves = wav)
geeglm_AR1_wav

# AR1
geeglm_AR1 <- geeglm(yvar ~ tvar, id = idvar, data = simdat, corstr = "ar1")
geeglm_AR1

# Independence
geeglm_ind <- geeglm(yvar ~ tvar, id = idvar, data = simdat, corstr = "independence")
geeglm_ind

# Exchangeable
geeglm_exch <- geeglm(yvar ~ tvar, id = idvar, data = simdat, corstr = "exchangeable")
geeglm_exch

# Unstructured
geeglm_unstr <- geeglm(yvar ~ tvar, id = idvar, data = simdat, corstr = "unstructured")
geeglm_unstr

# User-defined
zcor <- genZcor(clusz = table(simdat$idvar), 
                waves = simdat$timeorder, 
                corstrv = 4)

geeglm_user <- geeglm(yvar ~ tvar, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor, corstr = "userdefined")

# User-defined, Toeplitz
zcor_toep     <- matrix(NA, nrow(zcor), 4)
zcor_toep[,1] <- apply(zcor[,c(1, 5, 8,10)], 1, sum)
zcor_toep[,2] <- apply(zcor[,c(2, 6, 9)], 1, sum)
zcor_toep[,3] <- apply(zcor[,c(3, 7)], 1, sum)
zcor_toep[,4] <- zcor[,4]

geeglm_toep <- geeglm(yvar ~ tvar, 
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

geeglm_fix <- geeglm(yvar ~ tvar, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor_fix, corstr = "fixed")
geeglm_fix

test_that("bread works", {
  
  expect_true(check_bread(geeglm_AR1_wav, cluster = simdatPerm$idvar, y = simdatPerm$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_AR1, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_ind, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_exch, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_unstr, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_user, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_toep, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_fix, cluster = simdat$idvar, y = simdat$yvar))
  
})


test_that("vcovCR works for clustering variables higher than id variable.", {
  
  # create higher level
  pair_id <- rep(1:nlevels(as.factor(simdat$idvar)), each = 3, length.out = nlevels(as.factor(simdat$idvar)))[as.factor(simdat$idvar)] # factor cluster
  dat_higher <- cbind(simdat, pair_id)
  dat_scramble <- dat_higher[sample(nrow(dat_higher)),]
  
  # cluster at higher level
  
  expect_is(vcovCR(geeglm_AR1, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_exch, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_ind, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_unstr, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_user, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_toep, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  expect_is(vcovCR(geeglm_fix, type = "CR2", cluster = dat_higher$pair_id), "vcovCR")
  
  V_AR1 <- vcovCR(geeglm_AR1, type = "CR2", cluster =  dat_higher$pair_id)
  V_AR1_wav <- vcovCR(geeglm_AR1_wav, type = "CR2", cluster =  dat_higher$pair_id)
  V_exch <- vcovCR(geeglm_exch, type = "CR2", cluster =  dat_higher$pair_id)
  V_ind <- vcovCR(geeglm_ind, type = "CR2", cluster =  dat_higher$pair_id)
  V_unstr <- vcovCR(geeglm_unstr, type = "CR2", cluster =  dat_higher$pair_id)
  V_user <- vcovCR(geeglm_user, type = "CR2", cluster =  dat_higher$pair_id)
  V_toep <- vcovCR(geeglm_toep, type = "CR2", cluster =  dat_higher$pair_id)
  V_fix <- vcovCR(geeglm_fix, type = "CR2", cluster =  dat_higher$pair_id)
  expect_is(V_AR1, "vcovCR")
  expect_is(V_AR1_wav, "vcovCR")
  expect_is(V_exch, "vcovCR")
  expect_is(V_ind, "vcovCR")
  expect_is(V_unstr, "vcovCR")
  expect_is(V_user, "vcovCR")
  expect_is(V_toep, "vcovCR")
  expect_is(V_fix, "vcovCR")
  
  # check that result does not depend on sort-order
  expect_error(vcovCR(geeglm_AR1, type = "CR2", cluster = dat_scramble$pair_id)) # Warning
  V_scramble <- vcovCR(update(geeglm_AR1, data = dat_scramble), 
                       type = "CR2", cluster = dat_scramble$pair_id)
  expect_equal(diag(V), diag(V_scramble), tol = 10^-6)
})


test_that("vcovCR options work for CR2", {
  
  CR2_AR1_wav <- vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR2")
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR2"), CR2_AR1_wav) # calling bread from sandwich
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR2", 
                      inverse_var = TRUE), CR2_AR1_wav)
  expect_false(identical(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = simdatPerm$idvar, 
                                inverse_var = FALSE), CR2_AR1_wav))
    
  target <- targetVariance(geeglm_AR1_wav, cluster = simdatPerm$idvar) 
  expect_equal(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = simdatPerm$idvar, 
                        target = target, inverse_var = TRUE), CR2_AR1_wav, ignore_attr = TRUE)
  attr(CR2_AR1_wav, "inverse_var") <- FALSE
  expect_equal(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = simdatPerm$idvar, 
                        target = target, inverse_var = FALSE), CR2_AR1_wav, ignore_attr = TRUE)
    
})
  
  
test_that("vcovCR options work for CR4", {
  CR4_AR1_wav <- vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4")
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4"), CR4_AR1_wav)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = TRUE), CR4_AR1_wav)
  expect_false(identical(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = FALSE), CR4_AR1_wav))
  
  target <- targetVariance(geeglm_AR1_wav)
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
  expect_true(check_CR(geeglm_unstr, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_user, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_toep, cluster = simdat$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_fix, cluster = simdat$idvar, vcov = "CR2"))

  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_AR1, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_ind, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_exch, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_unstr, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_user, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_toep, cluster = simdat$idvar, vcov = "CR4"))
  expect_true(check_CR(geeglm_fix, cluster = simdat$idvar, vcov = "CR4"))
  
})



CR_types <- paste0("CR",0:4)

test_that("Order doesn't matter.", {
  
  check_sort_order(geeglm_ind, dat = simdat, cluster = simdat$idvar,
                   tol = 10^-4, tol2 = 10^-3, tol3 = 10^-3)
  
}) # Error
  
  
test_that("clubSandwich works with dropped observations", {
  dat_miss <- simdat
  dat_miss$yvar[sample.int(nrow(simdat), size = round(nrow(simdat) / 10))] <- NA
  dat_dropped <- geeglm(yvar ~ tvar, id = idvar, 
                       data = dat_miss, corstr = "independence")
  dat_complete <- geeglm(yvar ~ tvar, id = idvar, 
                                 data = simdat, corstr = "independence")
  
  CR_drop <- lapply(CR_types, function(x) vcovCR(dat_dropped, cluster = simdat$idvar, type = x))
  CR_complete <- lapply(CR_types, function(x) vcovCR(dat_complete, cluster = simdat$idvar, type = x))
  expect_equal(CR_drop, CR_complete)
  
  test_drop <- lapply(CR_types, function(x) coef_test(dat_dropped, cluster = simdat$idvar, vcov = x, test = "All", p_values = FALSE))
  test_complete <- lapply(CR_types, function(x) coef_test(dat_complete, cluster = simdat$idvar, vcov = x, test = "All", p_values = FALSE))
  expect_equal(test_drop, test_complete)
})


# add check.attributes = F to expect_equal
# vcovCR for geeglm is using bread function from sandwich package
# userdefined  genZcor function generate output with attributes at the end, nrow is not correct
# "bread works": not for the unstructured
# "vcovCR works for clustering variables higher than id variable." & "Order doesn't matter." has problem.
# "clubSandwich works with dropped observations"

