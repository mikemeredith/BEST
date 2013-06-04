BESTpower <-
function( BESTobj, N1, N2, credMass=0.95, ROPEm, ROPEsd, ROPEeff,
                     maxHDIWm, maxHDIWsd, maxHDIWeff, compValm=0, nRep=200,
                     mcmcLength=10000, saveName=NULL,
                     showFirstNrep=0 ) {
  # This function estimates power.
 
  # Sanity checks:
  if(!inherits(BESTobj, "data.frame"))
    stop("BESTobj is not a valid BEST object")
  if(ncol(BESTobj) == 3 && all(colnames(BESTobj) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
  } else if (ncol(BESTobj) == 5 && all(colnames(BESTobj) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("BESTobj is not a valid BEST object")
  }
  chainLength = NROW( BESTobj )
  if(chainLength < nRep)
    stop(paste("BESTobj does not have enough values; needs", nRep))
  if(credMass <= 0 || credMass >= 1)
    stop("credMass must lie between 0 and 1.")
  if(missing(N1))
    N1 <- length(attr(BESTobj, "data")$y1)
  N1 <- rep(N1, length.out=nRep)
  #if(!oneGrp && length(N2) < nRep)
  if(!oneGrp) {
    if(missing(N2))
      N2 <- length(attr(BESTobj, "data")$y2)
    N2 <- rep(N2, length.out=nRep)
  }


  # Deal with missing or invalid arguments for criteria:
  wanted <- rep(TRUE, 12)
  if(missing(ROPEm) || length(ROPEm) != 2)  {
    wanted[1:3] <- FALSE
    ROPEm <- c(NA_real_, NA_real_)
  }
  if(missing(ROPEsd) || length(ROPEsd) != 2)  {
    wanted[5:7] <- FALSE
    ROPEsd <- c(NA_real_, NA_real_)
  }
  if(missing(ROPEeff) || length(ROPEeff) != 2)  {
    wanted[9:11] <- FALSE
    ROPEeff <- c(NA_real_, NA_real_)
  }
  if(missing(maxHDIWm) || maxHDIWm <= 0)  {
    wanted[4] <- FALSE
    maxHDIWm <- NA_real_
  }
  if(missing(maxHDIWsd) || maxHDIWsd <= 0)  {
    wanted[8] <- FALSE
    maxHDIWsd <- NA_real_
  }
  if(missing(maxHDIWeff) || maxHDIWeff <= 0)  {
    wanted[12] <- FALSE
    maxHDIWeff <- NA_real_
  }
  if(!any(wanted))
    stop("No valid criteria set.")

  # Select thinned steps in chain for posterior predictions:
  stepIdxVec = seq( 1 , chainLength , floor(chainLength/nRep) )[1:nRep]
  paramDF <- BESTobj[stepIdxVec, ]
  
  goalTally <- numeric(12) 
  power <- matrix(NA, 12, 3)
  colnames(power) <- c("mean", "CrIlo", "CrIhi") # "CrI", cos too many HDIs already!
  rownames(power) <- c(
    "  mean:   HDI > ROPE",
    "  mean:   HDI < ROPE",
    "  mean:  HDI in ROPE",
    "  mean: HDI width ok",
    "    sd:   HDI > ROPE",
    "    sd:   HDI < ROPE",
    "    sd:  HDI in ROPE",
    "    sd: HDI width ok",
    "effect:   HDI > ROPE",
    "effect:   HDI < ROPE",
    "effect:  HDI in ROPE",
    "effect: HDI width ok")

  for (i in 1:nRep) {
    cat( "\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n" )
    cat( paste( "Power computation: Simulated Experiment" , i , "of" , 
                nRep , ":\n\n" ) )
    flush.console()
    # Get parameter values for this simulation:
    if(oneGrp) {
      mu1Val = paramDF[i,"mu"]
      sigma1Val = paramDF[i,"sigma"]
    } else {
      mu1Val = paramDF[i,"mu1"]
      mu2Val = paramDF[i,"mu2"]
      sigma1Val = paramDF[i,"sigma1"]
      sigma2Val = paramDF[i,"sigma2"]
    }
    nuVal = paramDF[i,"nu"]
    # Generate simulated data:
    y1 <- rt(N1[i], df=nuVal) * sigma1Val + mu1Val    
    y2 <- if(oneGrp) NULL else rt(N2[i], df=nuVal) * sigma2Val + mu2Val    
    # Get posterior for simulated data:
    simChain <- BESTmcmc( y1, y2, numSavedSteps=mcmcLength, thinSteps=1)
    if (i <= showFirstNrep ) { 
      x11()  # Changed 09-03-2013
      # dev.new() # Doesn't work in Rstudio ??
      plotAll(simChain, ROPEm=ROPEm, ROPEsd=ROPEsd, ROPEeff=ROPEeff,
              compValm=compValm) 
    }
    #simChain <- as.matrix(Bout)
    # Get the HDIs for each parameter:
    if(oneGrp) {
      HDIm <- hdi(simChain$mu, credMass=credMass) 
      HDIsd <- hdi(simChain$sigma, credMass=credMass) 
      mu0 <- if(is.null(compValm)) 0 else compValm
      HDIeff <- hdi((simChain$mu - mu0) / simChain$sigma,
                    credMass=credMass)
    } else {
      HDIm <- hdi(simChain$mu1 - simChain$mu2, credMass=credMass) 
      HDIsd = hdi(simChain$sigma1 - simChain$sigma2,
        credMass=credMass) 
      HDIeff = hdi(( simChain$mu1 - simChain$mu2 ) /
        sqrt( ( simChain$sigma1^2 + simChain$sigma2^2 ) / 2 ),
        credMass=credMass) 
    }
    # Assess which goals were achieved and tally them:
    goalTally <- goalTally + c(
      HDIm[1] > ROPEm[2],
      HDIm[2] < ROPEm[1],
      HDIm[1] > ROPEm[1] &  HDIm[2] < ROPEm[2],
      HDIm[2] - HDIm[1] < maxHDIWm,
      HDIsd[1] > ROPEsd[2],
      HDIsd[2] < ROPEsd[1],
      HDIsd[1] > ROPEsd[1] &  HDIsd[2] < ROPEsd[2],
      HDIsd[2] - HDIsd[1] < maxHDIWsd,
      HDIeff[1] > ROPEeff[2],
      HDIeff[2] < ROPEeff[1],
      HDIeff[1] > ROPEeff[1] &  HDIeff[2] < ROPEeff[2],
      HDIeff[2] - HDIeff[1] < maxHDIWeff )
      
    s1 = 1 + goalTally
    s2 = 1 + i - goalTally
    power[,1] = s1/(s1+s2)
    for ( j in which(wanted)) {
      power[j, 2:3] = hdi( qbeta , shape1=s1[j] , shape2=s2[j] )
    }
    cat( "\nAfter", i, "Simulated Experiments, Posterior Probability
       of meeting each criterion is (mean and 95% CrI):\n" )
    print(round(power[wanted, ], 3))
    flush.console()
    if(!is.null(saveName))
      save( i , power , file=saveName )
  }
  return(invisible(power))
}
