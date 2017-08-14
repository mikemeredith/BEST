# Post. pred. plots modified 14 Aug 2017 to show the data as an ordinary histogram.

plotAll <-
function(BESTobj, credMass=0.95,
                    ROPEm=NULL, ROPEsd=NULL, ROPEeff=NULL,
                    compValm=0, compValsd=NULL, compValeff=0,
                    showCurve=FALSE, ...) {
  # This function plots the posterior distribution (and data). It does not
  #   produce a scatterplot matrix; use pairs(...) for that.
  # Description of arguments:
  # BESTobj is BEST object of the type returned by function BESTmcmc.
  # ROPEm is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of means.
  # ROPEsd is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the difference of standard deviations.
  # ROPEeff is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE on the effect size.
  # showCurve if TRUE the posterior should be displayed as a fitted density curve
  #   instead of a histogram (default).

  # Sanity checks:
  if(!inherits(BESTobj, "data.frame"))
    stop("BESTobj is not a valid BEST object")
  if(ncol(BESTobj) == 3 && all(colnames(BESTobj) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
    colnames(BESTobj) <- c("mu1","nu","sigma1")
  } else if (ncol(BESTobj) == 5 && all(colnames(BESTobj) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("BESTobj is not a valid BEST object")
  }

  # Set up window and layout
  # ------------------------
  # windows(width=6.0,height=8.0) # Better to use default plot window.
  oldpar <- par(mar=c(3.5,3.5,2.5,0.5), mgp=c(2.25,0.7,0), "mfrow") 
    on.exit(par(oldpar))
  if(oneGrp) {
    layout( matrix( c(3,3,4,4,5,5, 1,1,1,1,2,2) , nrow=6, ncol=2 , byrow=FALSE ) )
  } else {
    layout( matrix( c(4,5,7,8,3,1,2,6,9,10) , nrow=5, byrow=FALSE ) )
  }

  # Do posterior predictive curve plots
  # -----------------------------------
  data <- attr(BESTobj, "data")

  # Select thinned steps in chain for plotting of posterior predictive curves:
  chainLength <- NROW( BESTobj )
  nCurvesToPlot <- 30
  stepIdxVec <- seq(1, chainLength, length.out=nCurvesToPlot)
  toPlot <- BESTobj[stepIdxVec, ]

  plotDataPPC(toPlot=toPlot, oneGrp=oneGrp, data=data)

  # Plot posterior distributions and their differences
  # --------------------------------------------------
  # Plot posterior distribution of parameter nu:
  plotPost( log10(BESTobj$nu) , col="skyblue" ,  credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote("log10("*nu*")") , cex.lab = 1.75 , showMode=TRUE ,
                  main="Normality" ) #  (<0.7 suggests kurtosis)

  # Plot posterior distribution of parameters mu1, mu2, and their difference:
  xlim <- range(BESTobj$mu1, BESTobj$mu2)
  if(oneGrp) {
    plotPost( BESTobj$mu1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve , ROPE=ROPEm, compVal=compValm,
                  xlab=bquote(mu) , main="Mean" , 
                  col="skyblue" )
  } else {
    plotPost( BESTobj$mu1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote(mu[1]) , main=paste("Group",1,"Mean") , 
                  col="skyblue" )
    plotPost( BESTobj$mu2 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                         showCurve=showCurve ,
                    xlab=bquote(mu[2]) , main=paste("Group",2,"Mean") , 
                    col="skyblue" )
    plotPost( BESTobj$mu1-BESTobj$mu2 , compVal=0 ,  showCurve=showCurve , credMass=credMass,
                    xlab=bquote(mu[1] - mu[2]) , cex.lab = 1.75 , ROPE=ROPEm ,
                    main="Difference of Means" , col="skyblue" )
  }

  # Plot posterior distribution of param's sigma1, sigma2, and their difference:
  xlim <- range(BESTobj$sigma1, BESTobj$sigma2)
  if(oneGrp) {
    plotPost(BESTobj$sigma1, xlim=xlim, cex.lab = 1.75, credMass=credMass,
                       showCurve=showCurve, ROPE=ROPEsd, compVal=compValsd,
                  xlab=bquote(sigma) , main="Std. Dev." , 
                  col="skyblue" , showMode=TRUE )
  } else {
    plotPost( BESTobj$sigma1 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                       showCurve=showCurve ,
                  xlab=bquote(sigma[1]) , main=paste("Group",1,"Std. Dev.") , 
                  col="skyblue" , showMode=TRUE )
    plotPost( BESTobj$sigma2 ,  xlim=xlim , cex.lab = 1.75 , credMass=credMass,
                         showCurve=showCurve ,
                    xlab=bquote(sigma[2]) , main=paste("Group",2,"Std. Dev.") , 
                    col="skyblue" , showMode=TRUE )
    plotPost( BESTobj$sigma1-BESTobj$sigma2 ,  credMass=credMass,
                         compVal=compValsd ,  showCurve=showCurve ,
                         xlab=bquote(sigma[1] - sigma[2]) , cex.lab = 1.75 , 
                         ROPE=ROPEsd ,
                 main="Difference of Std. Dev.s" , col="skyblue" , showMode=TRUE )
  }

  # Plot effect size
  # ----------------
  # Effect size for 1 group:
  if(oneGrp)  {
  effectSize = ( BESTobj$mu1 - compValm ) / BESTobj$sigma1
  plotPost( effectSize , compVal=compValeff ,  ROPE=ROPEeff , credMass=credMass,
                 showCurve=showCurve ,
                 xlab=bquote( (mu-.(compValm)) / sigma ),
                 showMode=TRUE , cex.lab=1.75 , main="Effect Size" , col="skyblue" )
  } else {
    # Plot of estimated effect size. Effect size is d-sub-a from 
    # Macmillan & Creelman, 1991; Simpson & Fitter, 1973; Swets, 1986a, 1986b.
    effectSize <- ( BESTobj$mu1 - BESTobj$mu2 ) / sqrt( ( BESTobj$sigma1^2 + BESTobj$sigma2^2 ) / 2 )
    plotPost( effectSize , compVal=compValeff ,  ROPE=ROPEeff , credMass=credMass,
                          showCurve=showCurve ,
                    xlab=bquote( (mu[1]-mu[2])
                      /sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
                showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  }
  # Or use sample-size weighted version:
  # Hedges 1981; Wetzels, Raaijmakers, Jakab & Wagenmakers 2009.
  # N1 = length(y1)
  # N2 = length(y2)
  # effectSize = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) )
  #                                    / (N1+N2-2) )
  # Be sure also to change BESTsummary function, above.
  # histInfo = plotPost( effectSize , compVal=0 ,  ROPE=ROPEeff ,
  #          showCurve=showCurve ,
  #          xlab=bquote( (mu[1]-mu[2])
  #          /sqrt((sigma[1]^2 *(N[1]-1)+sigma[2]^2 *(N[2]-1))/(N[1]+N[2]-2)) ),
  #          showMode=TRUE , cex.lab=1.0 , main="Effect Size" , col="skyblue" )
  return(invisible(NULL))
}
