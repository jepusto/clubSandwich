
# This script has all the support code to run a large simulation to look at
# stability of the different SE estimators and true performances of the ATE
# estimators.
#
# The simulation script runs a bunch of scenarios and saves the massive amount
# of results to a file to be processed.
#
# This script has the code to run a single trial, etc.

library( dplyr )
library( blkvar )
library( purrr )
library( multiwayvcov )

if ( FALSE ) {
  # for testing crashing in the simulator
  my_compare_methods =function(... ){
    if ( sample(5,1) == 1 ) {
      asdfadfga
    } else {
      compare_methods( ... )
    }
  }
}


scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}


get.estimates = function( df, ... ) {
  safe_comp = safely( compare_methods )
  ests = safe_comp( Yobs, Z, sid, data=df, ... )
  if ( is.null( ests[[2]] ) ) {
    ests = ests[[1]]
    ests = rename( ests, ATE.hat = tau,
                   SE.hat = SE )
    
  } else {
    cat( "\nCaught error:\n" )
    print( ests[[2]] )
    ests = data.frame( ATE.hat = c(NA),
                       SE.hat = c(NA) )    
  }
  ests
}

# Reassign treatment to do finite sample inference.
rerandomize.data = function( dat ) {
  dat = dat %>% group_by( sid ) %>%
    mutate( Z = sample( Z ) ) %>% ungroup()
  dat = mutate( dat, Yobs = ifelse( Z, Y1, Y0 ) )
  dat
}





describe.data = function( df ) {
  
  tau.S = attr( df, "tau.S")
  
  df = rerandomize.data( df )
  
  params = df %>% group_by( sid ) %>%
    summarise( ATE = mean( Y1 ) - mean( Y0 ),
               n = n(),
               ATE.hat = mean( Y1[Z==1] ) - mean( Y0[Z==0] ),
               p = mean( Z == 1 ) )
  
  ATE.person = with( df, mean( Y1 - Y0 ) )
  
  sd.p = sd( params$p )
  sd.n = sd( params$n )
  
  if (FALSE) {
    M = lm( Yobs ~ Z * sid - 0 - Z, data=df )
    summary( M )
    (summary(M))$r.squared
    sd( params$ATE - params$ATE.hat )
    sd( params$ATE )
    var( params$ATE )
    var( params$ATE - params$ATE.hat )
    sd( df$Y0 )
  }
  
  sstat = params %>% summarise( ATE.site = mean( ATE ),
                                sd.ATE = sd( ATE ),
                                sd.ATE.hat = sd( ATE.hat ),
                                sd.n = sd( n ),
                                sd.p = sd( p ),
                                cor.p.Bhat = ifelse( sd.p > 0, cor( ATE.hat, p ), NA ),
                                cor.n.Bhat = ifelse( sd.n > 0, cor( ATE.hat, n ), NA ),
                                cor.p = ifelse( sd.p > 0, cor( ATE, p ), NA ),
                                cor.n = ifelse( sd.n > 0, cor( ATE, n ), NA ) )
                  
  sstat = mutate( sstat,
                  ATE.finite.person=ATE.person )
  sstat
}


# Generate a multisite trial, estimate ATE using all our methods, and then also
# calculate the true ATE (person and site weighted).
# @return Dataframe with the estimates along with the true baseline values.
single.MLM.trial = function( n.bar, J, tau.11.star, dependence = FALSE, proptx.dependence = FALSE, variable.n = FALSE, variable.p = FALSE,
                             ATE.superpop = 0.2, ICC = 0.20, p.tx = 0.65, n.runs = 3,
                             just.describe.data = FALSE, ... ) {
  
  df = gen.dat.no.cov( n.bar=n.bar, 
                       J=J,
                       tau.11.star = tau.11.star,
                       gamma.10 = ATE.superpop,
                       ICC = ICC,
                       p = p.tx,
                       variable.n = variable.n,
                       variable.p = variable.p,
                       finite.model = FALSE,
                       size.impact.correlate = dependence,
                       proptx.impact.correlate = proptx.dependence,
                       correlate.strength = 0.50,
                       size.ratio = 0.60
                      )
  
  if ( just.describe.data ) {
    return( describe.data( df ) )
  }
  
  tau.S = attr( df, "tau.S")
  
  params = df %>% group_by( sid ) %>%
    summarise( ATE = mean( Y1 ) - mean( Y0 ),
               n = n() )
  ATE.person = with( df, mean( Y1 - Y0 ) )
  ATE.site = mean( params$ATE )
  
  ests = plyr::rdply( n.runs, {
    df = rerandomize.data( df )
    get.estimates( df, ... ) 
  }, .id="subrun" )
  
  ests = mutate( ests,
                 ATE.finite.person=ATE.person, 
                 ATE.finite.site=ATE.site )
  ests    
}



# Testing code
if ( FALSE ) {
  df = gen.dat.no.cov( n.bar=200, J=30,
                       tau.11.star = 0.1^2,
                       ICC = 0.20,
                       p = 0.70,
                       variable.n = TRUE,
                       variable.p = TRUE,
                       size.impact.correlate = TRUE,
                       proptx.impact.correlate = TRUE,
                       finite.model = FALSE )
  tau.S = attr( df, "tau.S")
  tau.S    
  
  sites = df %>% group_by( sid ) %>% 
    summarise( n = n(),
               p.Z = mean( Z ),
               ATE.hat = mean( Yobs[Z==1] ) - mean( Yobs[Z==0] ) )
  
  head( sites )
  qplot( sites$n )
  summary( sites$n )
  qplot( sites$p.Z )
  summary( sites$p.Z )
  qplot( sites$n, sites$ATE.hat )
  cor( sites$n, sites$ATE.hat )
  
  qplot( sites$p.Z, sites$ATE.hat )
  cor( sites$p.Z, sites$ATE.hat )
  
  rst = single.MLM.trial( 20, 10, 0.2^2, dependence = TRUE, proptx.dependence = TRUE )
  head( rst )
  table( rst$subrun )
}

#single.MLM.trial( 20, 10, 0.2^2 )



# Run a simulation with R trials and return all the simulation runs as a large dataframe.
run.scenario = function( J, n.bar, tau, dependence, proptx.dependence, variable.n, variable.p, ATE = 0.20, ICC = 0.20, p.tx = 0.65, R = 10, n.runs=3, 
                         ID=NULL, .progress="none", ... ) {
  ptm = proc.time()
  
  #  n = n.bar * J
  #  e.time = 1.13^(J/5) + 1.0183^(n/100)
  
  #  R.adj = round( R * pmin( 1, 3/e.time ) )
  #if ( J > 40 ) {
  #  R.adj = R / 2
  #} else {
  #  R.adj = R
  #}
  #R.adj = R
  
  if ( is.null( ID ) ) {
    ID = -999
  }
  
  scat( "Running scenario %d\n\tJ=%d\tn.bar=%d\ttau=%.2f\tATE=%.2f\tICC=%.2f\tprop tx=%.2f\tdependence=%s\tprop dep=%s\tR=%d\n", 
        ID, J, n.bar,tau, ATE, ICC, p.tx, dependence, proptx.dependence, R)
  
  rps = plyr::rdply( R,  single.MLM.trial( n.bar=n.bar, J=J, tau.11.star=tau, 
                                               dependence = dependence, proptx.dependence = proptx.dependence,
                                               variable.n = variable.n, variable.p = variable.p,
                                               ATE.superpop=ATE, ICC=ICC, p.tx=p.tx, n.runs=n.runs, ... ),
                     .id="run", .progress="text" )
  
  rps$subrun = paste0( rps$run, ".", rps$subrun )
  
  scat("**\n**\tTotal time elapsed:\n")
  tot.time = proc.time() - ptm
  print(tot.time)
  sim.per.min = R / (tot.time["elapsed"] / 60)
  scat("Simulations per minute = %.2f\n", sim.per.min )
  rps$sim.per.min = sim.per.min
  rps
}


# Testing run.scenario()
if ( FALSE ) {
  rr = run.scenario( J = 4, n.bar = 16, tau = 0.3^2, dependence=TRUE, proptx.dependence = TRUE, 
                     variable.n = TRUE, variable.p = TRUE,
                     ATE = 0.2, R = 10 )
  head( rr )  
  table( rr$run )
  table(  table( rr$run ) )
  table(  table( rr$subrun ) )
  tt = table( rr$subrun )
  tt
  tt[ tt < 10 ]
  filter( rr, subrun %in% c( "1.2", "4.1", "4.2", "2.1" ) )
  
}



#### Make scenario configurations #####

make.scenario.list = function( group = "main" ) {
  
  if ( group == "main" ) {
    scenarios = expand.grid( J = c( 80, 40, 20, 10 ),
                             n.bar = c( 8000, 4000, 2000, 1000 ), # put in totals here
                             dependence = c( 1, 0 ),
                             proptx.dependence = c( 1, 0, -1 ), 
                             variable.n = c( TRUE, FALSE ),
                             variable.p = c( TRUE, FALSE ),
                             ATE = 0.20,
                             tau = c( 0, 0.1, 0.20 )^2,
                             p.tx = 0.65,
                             ICC = 0.20 )
  } else if ( simstudy == "large" ) {
    scenarios = expand.grid( J = c( 80, 40, 20 ),
                             n.bar = c( 8000, 4000, 2000 ), # put in totals here
                             dependence = c( 1, 0 ),
                             proptx.dependence = c( 1, 0, -1 ), 
                             variable.n = c( TRUE, FALSE ),
                             variable.p = c( TRUE, FALSE ),
                             ATE = 0.20,
                             tau = c( 0, 0.1, 0.20 )^2,
                             p.tx = 0.70,
                             ICC = 0.60 )
  }
  
  nrow( scenarios )
  
  # drop redundant scenarios (dependence irrelevant if quantaty not varying)
  scenarios = filter( scenarios, (variable.n == TRUE) | (dependence == 0) )
  scenarios = filter( scenarios, (variable.p == TRUE) | (proptx.dependence == 0) )
  nrow( scenarios )
  
  # drop redundant/not well defined scenarios (dependence irrelevant if no cross site variation )
  scenarios = filter( scenarios, !( (tau == 0) & (dependence != 0) ) )
  scenarios = filter( scenarios, !( (tau == 0) & (proptx.dependence != 0) ) )
  nrow( scenarios )
  
  table( scenarios$variable.n, scenarios$dependence )
  table( scenarios$variable.p, scenarios$proptx.dependence )
  
  table( scenarios$variable.n, scenarios$tau )
  
  # convert n.bar from totals to people per site
  scenarios = mutate( scenarios, n.bar = round( n.bar / J ) )
  scenarios = as.tibble( scenarios )
  scenarios
  nrow( scenarios )
  
  scenarios$ID = 1:nrow( scenarios )
  
  scenarios  
}


get.scenario.by.id = function( index, group = "main" ) {
  scenarios = make.scenario.list( group=group )
  scenarios[ index, ]
}



if ( FALSE ) {
  sc = make.scenario.list()
  nrow( sc )
  head( sc )
  table( sc$variable.n, sc$variable.p )
}


