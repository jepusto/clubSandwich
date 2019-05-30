##
## Code from tracking down bugs with clubsandwich and cluster robust SEs
##

library( tidyverse )
library( blkvar )

# Mini sim study to look at club sandwich and the cluster robust SE stuff
source("auxilliary/simulation_support_code.R" )

set.seed(20190530)

lmsim = run.scenario(
  J = 10, 
  n.bar = 100, 
  tau=0.2^2, 
  dependence = 0, 
  proptx.dependence = 0, 
  variable.n = TRUE, 
  variable.p = FALSE, 
  include.MLM=FALSE, 
  include.block=FALSE, 
  n.runs=1, 
  R = 1000, 
  .progress="text")

table( lmsim$method )

clubsim = filter( lmsim, method %in% c( "FE-CR", "FE-Club" ) )
head( clubsim )

# Look at club sandwich
rst = clubsim %>% 
  group_by( method ) %>% 
  summarise(
    SE = sd( ATE.hat ),
    V = var(ATE.hat),
    E.SE.hat = mean( SE.hat ),
    E.V.hat = mean(SE.hat^2),
    med.SE.hat = median(SE.hat),
    sd.SE.hat = sd(SE.hat),
    rat_SE = E.SE.hat / SE,
    rat_V = E.V.hat / V
  )
rst




# Quick simulation check:
#
# Does hand-rolled code give same results as the clubsandwich package code?
# Aside: major time difference in running!

library(clubSandwich)
library(plm)

check = replicate( 10, {
  df = gen.dat.no.cov( n.bar=400, J=20,
                       tau.11.star = 0,
                       ICC = 0.20,
                       p = 0.70,
                       variable.n = TRUE,
                       variable.p = TRUE,
                       size.impact.correlate = TRUE,
                       proptx.impact.correlate = TRUE,
                       finite.model = FALSE )
  head( df )
  
  cm = compare_methods( Yobs, Z, sid, data=df, include.MLM = FALSE, include.block = FALSE, include.DB = FALSE )
  cm
  
  M0 = lm( Yobs ~ 0 + Z + sid, data=df )
  system.time(ct <- coef_test( M0, vcov="CR2", cluster = df$sid, coefs=c("Z" ) ))
  ct
  
  PLM0 = plm(Yobs ~ Z, data = df, effect = "individual", model = "within", index = "sid")
  system.time(plm_ct <- coef_test(PLM0, vcov = "CR2", ignore_FE = TRUE))
  plm_ct
  
  SEs = c( hand = filter( cm, method=="FE-Club" )$SE,
           package_lm = ct$SE,
           package_plm = plm_ct$SE)
  if ( diff( range( SEs ) ) > 0.05 ) {
    browser()
  }
  SEs
} )

check 
check[1,] - check[2,]


