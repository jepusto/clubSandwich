
# skip("Just for debugging purposes.")

library(geepack)

data(dietox)
dietox$Cu <- as.factor(dietox$Cu)

cor_fix <- 0.4^as.matrix(dist(1:12))

zcor_fix <- fixed2Zcor(cor_fix, id = dietox$Pig, waves = dietox$Time)

suppressWarnings(
  mod <- geeglm(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)), 
                data=dietox, 
                id=Pig, 
                family=poisson("identity"), 
                zcor = zcor_fix, corstr = "fixed")
)

get_correlations <- function(obj) {
  i <- 1L
  while(!identical(parent.frame(i), .GlobalEnv)) {
    cat("\n i:", i, "||", ls(envir = parent.frame(n = i)))
    i <- i + 1L
  }
  cat("\n i:", i, "||", ls(envir = parent.frame(n = i)))
  
  cat("\n", "Envir: ", find(as.character(obj$call$zcor)), find(as.character(obj$call$zcor), numeric = TRUE), "\n")
  eval(obj$call$zcor, envir = parent.frame())
  
}

something <- function(obj) {
  something <- get_correlations(obj)
  dim(something)
}

get_correlations(mod)
something(mod)
# v <- targetVariance(mod, cluster = dietox$Pig)
# w <- weightMatrix(mod, cluster = dietox$Pig)

# test_that("weightMatrix() works for geeglm with corstr = 'fixed'.", {
# 
#   corr_vec <- get_correlations(mod)
#   Vmat <- targetVariance(mod, cluster = dietox$Pig)
#   Wmat <- weightMatrix(mod, cluster = dietox$Pig)
#   expect_length(corr_vec, length(zcor_fix))
# 
# })
