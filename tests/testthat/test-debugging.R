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
  for (i in 1:10) cat("\n", "i:", i, "||", ls(envir = parent.frame(n = i)))
  cat("\n", "Envir: ", find(as.character(obj$call$zcor)), find(as.character(obj$call$zcor), numeric = TRUE), "\n")
  zcor <- eval(obj$call$zcor, enclos = parent.frame())

  return(zcor)  
}

test_that("targetVariance works for geeglm with corstr = 'fixed'.", {

  Vmat <- targetVariance(mod, cluster = dietox$Pig)
  Wmat <- weightMatrix(mod, cluster = dietox$Pig)
  corr_vec <- get_correlations(mod)
  expect_length(corr_vec, length(zcor_fix))
  
})
