plotDataPPC <-
function(y, mu, sigma,  nu, xVec, stepIdxVec, maxY) {
  # Does the plots of posterior predictive curves; called by plotAll;
  #  not exported.
  # Does not do title or sample size: those are added later.
  stepIdx <- 1
  plot( xVec , dt( (xVec-mu[stepIdxVec[stepIdx]])/sigma[stepIdxVec[stepIdx]] , 
                   df=nu[stepIdxVec[stepIdx]] )/sigma[stepIdxVec[stepIdx]] , 
        ylim=c(0,maxY) , cex.lab=1.75 ,
        type="l" , col="skyblue" , lwd=1 , xlab="y" , ylab="p(y)" )
  for ( stepIdx in 2:length(stepIdxVec) ) {
    lines(xVec, dt( (xVec-mu[stepIdxVec[stepIdx]])/sigma[stepIdxVec[stepIdx]] , 
                      df=nu[stepIdxVec[stepIdx]] )/sigma[stepIdxVec[stepIdx]] , 
           type="l" , col="skyblue" , lwd=1 )
  }
  histBinWd <- median(sigma)/2
  histCenter <- mean(mu)
  histBreaks <- sort( c( seq( histCenter-histBinWd/2 , min(xVec)-histBinWd/2 ,
                             -histBinWd ),
                        seq( histCenter+histBinWd/2 , max(xVec)+histBinWd/2 ,
                             histBinWd ) ) )
  if(!is.null(y)) {
    histInfo <- hist( y, plot=FALSE , breaks=histBreaks )
    PlotMat <- cbind(histInfo$mids, histInfo$density)
    PlotMat[histInfo$density == 0] <- NA
    points( PlotMat, type="h" , lwd=3 , col="red" )
  }
}
