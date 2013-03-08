BESTmcmc <-
function( y1, y2=NULL,
    numSavedSteps=1e5, thinSteps=1, burnInSteps = 1000) { 
  # This function generates an MCMC sample from the posterior distribution.
  # y1, y2 the data vectors; y2=NULL if only one group.
  # Returns an mcmc.list object, not a matrix, with class 'BEST'.
  
  require(rjags)

  #------------------------------------------------------------------------------
  modelFile <- file.path(tempdir(), "BESTmodel.txt")
  # THE MODEL.
  if(is.null(y2)) {
    modelString = "
    model {
      for ( i in 1:Ntotal ) {
        y[i] ~ dt( mu , tau , nu )
      }
      mu ~ dnorm( muM , muP )
      tau <- 1/pow( sigma , 2 )
      sigma ~ dunif( sigmaLow , sigmaHigh )
      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }
    " # close quote for modelString
  } else {
    modelString <- "
    model {
      for ( i in 1:Ntotal ) {
        y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )
      }
      for ( j in 1:2 ) {
        mu[j] ~ dnorm( muM , muP )
        tau[j] <- 1/pow( sigma[j] , 2 )
        sigma[j] ~ dunif( sigmaLow , sigmaHigh )
      }
      nu <- nuMinusOne+1
      nuMinusOne ~ dexp(1/29)
    }
    " # close quote for modelString
  }   
  # Write out modelString to a text file
  writeLines( modelString , con=modelFile )
  
  #------------------------------------------------------------------------------
  # THE DATA.
  # Load the data:
  y = c( y1 , y2 ) # combine data into one vector
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    Ntotal = Ntotal ,
    muM = mean(y) ,
    muP = 0.000001 * 1/sd(y)^2 ,
    sigmaLow = sd(y) / 1000 ,
    sigmaHigh = sd(y) * 1000 
  )
  if(!is.null(y2)) # create group membership code
    dataList$x <- c( rep(1,length(y1)) , rep(2,length(y2)) )

  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  if(is.null(y2)) {
    mu = mean(y1)
    sigma = sd(y1)
  } else {
    mu = c( mean(y1) , mean(y2) )
    sigma = c( sd(y1) , sd(y2) )
  }
  # Regarding initial values in next line: (1) sigma will tend to be too big if 
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  initsList = list( mu = mu , sigma = sigma , nuMinusOne = 4 )
  
  #------------------------------------------------------------------------------
  # RUN THE CHAINS

  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  nChains = 3 
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  cat( "Setting up the JAGS model...\n" ) ; flush.console()
  jagsModel = jags.model( modelFile, data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" ) ; flush.console()
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" ) ; flush.console()
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  #------------------------------------------------------------------------------
 
  # mcmcChain = as.matrix( codaSamples )
  class(codaSamples) <- c("BEST", class(codaSamples))
  return( codaSamples )
}
