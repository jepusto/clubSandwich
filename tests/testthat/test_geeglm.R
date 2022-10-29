context("geeglm objects")
set.seed(202201028)

library(geepack)

# Simulating a dataset
timeorder <- rep(1:5, 10)
tvar      <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:10, each=5)
uuu   <- rep(rnorm(10), each=5)
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
nrow(zcor)
clusz <- table(simdat$idvar)
sum(clusz * (clusz - 1) / 2)

geeglm_user <- geeglm(yvar ~ tvar, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor, corstr = "userdefined")

zcor_toep     <- matrix(NA, nrow(zcor), 4)
zcor_toep[,1] <- apply(zcor[,c(1, 5, 8,10)], 1, sum)
zcor_toep[,2] <- apply(zcor[,c(2, 6, 9)], 1, sum)
zcor_toep[,3] <- apply(zcor[,c(3, 7)], 1, sum)
zcor_toep[,4] <- zcor[,4]

geeglm_toep <- geeglm(yvar ~ tvar, 
                      id = idvar, waves = timeorder, 
                      data = simdat, 
                      zcor = zcor_toep, corstr = "userdefined")



test_that("bread works", {
  
  expect_true(check_bread(geeglm_AR1_wav, cluster = simdatPerm$idvar, y = simdatPerm$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_AR1, cluster = simdat$idvar, y = simdat$yvar, tol = 1e-5))
  expect_true(check_bread(geeglm_ind, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_exch, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_unstr, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_user, cluster = simdat$idvar, y = simdat$yvar))
  expect_true(check_bread(geeglm_toep, cluster = simdat$idvar, y = simdat$yvar))
  
})


test_that("vcovCR works for clustering variables higher than id variable.", {
  
  # create higher level
  pair_id <- rep(1:nlevels(as.factor(simdat$idvar)), each = 3, length.out = nlevels(as.factor(simdat$idvar)))[as.factor(simdat$idvar)] # factor cluster
  
  re_order <- sample(nrow(simdat))
  dat_scramble <- simdat[re_order,]
  pair_scramble <- pair_id[re_order]
  
  # cluster at higher level
  expect_is(vcovCR(geeglm_ind, type = "CR2", cluster = pair_id), "vcovCR")
  V <- vcovCR(geeglm_ind, type = "CR2", cluster = pair_id)
  expect_is(V, "vcovCR")
  
  expect_error(vcovCR(geeglm_ind, type = "CR2", cluster = pair_scramble)) # No error
  
  # check that result does not depend on sort-order
  V_scramble <- vcovCR(update(geeglm_ind, data = dat_scramble), 
                       type = "CR2", cluster = pair_scramble)
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
                        target = target, inverse_var = TRUE), CR2_AR1_wav, check.attributes = FALSE)
  attr(CR2_AR1_wav, "inverse_var") <- FALSE
  expect_equal(vcovCR(geeglm_AR1_wav, type = "CR2", cluster = simdatPerm$idvar, 
                        target = target, inverse_var = FALSE), CR2_AR1_wav, check.attributes = FALSE)
    
})
  
  
test_that("vcovCR options work for CR4", {
  CR4_AR1_wav <- vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4")
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4"), CR4_AR1_wav)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = TRUE), CR4_AR1_wav)
  expect_false(identical(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", inverse_var = FALSE), CR4_AR1_wav))
  
  target <- targetVariance(geeglm_AR1_wav)
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", 
                      target = target, inverse_var = TRUE), CR4_AR1_wav, check.attributes = FALSE)
  attr(CR4_AR1_wav, "inverse_var") <- FALSE
  expect_equal(vcovCR(geeglm_AR1_wav, cluster = simdatPerm$idvar, type = "CR4", 
                      target = target, inverse_var = FALSE), CR4_AR1_wav, check.attributes = FALSE)
})


test_that("CR2 and CR4 are target-unbiased", {
  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar, vcov = "CR2"))
  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar,vcov = "CR2"))
  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar,vcov = "CR4"))
  expect_true(check_CR(geeglm_AR1_wav, cluster = simdatPerm$idvar,vcov = "CR4"))
})

test_that("get_data works.", {
  re_order <- sample(nrow(simdatPerm))
  sim_scramble <- simdatPerm[re_order,]
  geeglm_AR1_wav_scramble <- geeglm(yvar ~ tvar, id = idvar, 
                        data = sim_scramble, 
                        corstr = "ar1", waves = wav)
  scramble_dat <- get_data(geeglm_AR1_wav_scramble)
  expect_equal(sim_scramble, scramble_dat)
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

