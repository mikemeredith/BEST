makeData <-
function( mu1 , sd1 , mu2=NULL , sd2=NULL , nPerGrp , 
                     pcntOut=0 , sdOutMult=2.0 ,
                     rnd.seed=NULL, showPlot=TRUE ) {
  # Auxilliary function for generating random values from a 
  # mixture of normal distibutions.
  oneGrp <- is.null(mu2) || is.null(sd2)
  if(!is.null(rnd.seed)){set.seed(rnd.seed)} # Set seed for random values.
  nOut = ceiling(nPerGrp*pcntOut/100)        # Number of outliers.
  nIn = nPerGrp - nOut                       # Number from main distribution.
  #
  y1 = rnorm(n=nIn)                     # Random values in main distribution.
  y1 = ((y1-mean(y1))/sd(y1))*sd1 + mu1 # Translate so exactly realize desired.
  y1Out = rnorm(n=nOut)                 # Random values for outliers
  y1Out=((y1Out-mean(y1Out))/sd(y1Out))*(sd1*sdOutMult)+mu1 # Realize exactly.
  y1 = c(y1,y1Out)                      # Concatenate main with outliers.
  #
  if(oneGrp) {
    y2 <- NULL
  } else {
    y2 = rnorm(n=nIn)                     # Repeat for second group.
    y2 = ((y2-mean(y2))/sd(y2))*sd2 + mu2
    y2Out = rnorm(n=nOut) 
    y2Out=((y2Out-mean(y2Out))/sd(y2Out))*(sd2*sdOutMult)+mu2
    y2 = c(y2,y2Out)
  }
  #
  if(showPlot) {
    # Set up window and layout:
    opar <- par(mfrow=c(1,1)) ; on.exit(par(opar))
    if(!oneGrp)  
      par(mfrow=2:1)  
      # layout(matrix(1:2,nrow=2))  
    histInfo = hist( y1 , main="Simulated Data" , col="pink2" , border="white" , 
                     xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
    text( max(c(y1,y2)) , max(histInfo$density) , 
          bquote(N==.(nPerGrp)) , adj=c(1,1) )
    xlow=min(histInfo$breaks)
    xhi=max(histInfo$breaks)
    xcomb=seq(xlow,xhi,length=1001)
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) +
      dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) , lwd=3 )
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1)*nIn/(nIn+nOut) , 
           lty="dashed" , col="blue", lwd=3)
    lines( xcomb , dnorm(xcomb,mean=mu1,sd=sd1*sdOutMult)*nOut/(nIn+nOut) , 
           lty="dashed" , col="red", lwd=3)
    if(!oneGrp)  {
      histInfo = hist( y2 , main="" , col="pink2" , border="white" , 
                       xlim=range(c(y1,y2)) , breaks=30 , prob=TRUE )
      text( max(c(y1,y2)) , max(histInfo$density) , 
            bquote(N==.(nPerGrp)) , adj=c(1,1) )
      xlow=min(histInfo$breaks)
      xhi=max(histInfo$breaks)
      xcomb=seq(xlow,xhi,length=1001)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) +
        dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) , lwd=3)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2)*nIn/(nIn+nOut) , 
             lty="dashed" , col="blue", lwd=3)
      lines( xcomb , dnorm(xcomb,mean=mu2,sd=sd2*sdOutMult)*nOut/(nIn+nOut) , 
             lty="dashed" , col="red", lwd=3)
    }
  }
  #
  return( list( y1=y1 , y2=y2 ) )
}
