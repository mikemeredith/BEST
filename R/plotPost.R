plotPost <-
function( paramSampleVec, credMass=0.95, compVal=NULL, ROPE=NULL, 
           HDItextPlace=0.7, yaxt="n", ylab="",
           xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
           xlim=range(compVal, paramSampleVec), 
           col="skyblue", border="white", showMode=FALSE,
           showCurve=FALSE, breaks=NULL, ... ) {

  # Does a plot for a single parameter. Called by plot.BEST but also exported.
  # Returns invisible(NULL).

  oldpar <- par(xpd=NA) ; on.exit(par(oldpar))

  # Plot histogram.
  if ( is.null(breaks) ) {
    by <- diff(hdi(paramSampleVec))/18
    breaks <- c( seq( from=min(paramSampleVec), to=max(paramSampleVec),
                   by=by), max(paramSampleVec) )
  }
  if ( !showCurve ) {
    histinfo <- hist( paramSampleVec, xlab=xlab, yaxt=yaxt, ylab=ylab,
                     freq=FALSE, border=border, col=col,
                     xlim=xlim, main=main, cex=cex, cex.lab=cex.lab,
                     breaks=breaks, ... )
  }
  if ( showCurve ) {
    histinfo <- hist( paramSampleVec, plot=FALSE )
    densCurve <- density( paramSampleVec, adjust=2 )
    plot( densCurve$x, densCurve$y, type="l", lwd=5, col=col, bty="n",
          xlim=xlim, xlab=xlab, yaxt=yaxt, ylab=ylab,
          main=main, cex=cex, cex.lab=cex.lab, ... )
  }
  cenTendHt <- 0.9*max(histinfo$density)
  cvHt <- 0.7*max(histinfo$density)
  ROPEtextHt <- 0.55*max(histinfo$density)
  # Display mean or mode:
  if ( showMode==FALSE ) {
      meanParam <- mean( paramSampleVec )
      text( meanParam, cenTendHt,
            bquote(mean==.(signif(meanParam,3))), adj=c(.5,0), cex=cex )
  } else {
      dres <- density( paramSampleVec )
      modeParam <- dres$x[which.max(dres$y)]
      text( modeParam, cenTendHt,
            bquote(mode==.(signif(modeParam,3))), adj=c(.5,0), cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    cvCol <- "darkgreen"
    pcgtCompVal <- round( 100 * sum( paramSampleVec > compVal )
                          / length( paramSampleVec ) , 1 )
     pcltCompVal <- 100 - pcgtCompVal
     lines( c(compVal,compVal), c(0.96*cvHt,0),
            lty="dotted", lwd=1, col=cvCol )
     text( compVal, cvHt,
           bquote( .(pcltCompVal)*"% < " *
                   .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ),
           adj=c(pcltCompVal/100,0), cex=0.8*cex, col=cvCol )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ropeCol <- "darkred"
     pcInROPE <- ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                          / length( paramSampleVec ) )
     lines( c(ROPE[1],ROPE[1]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol )
     lines( c(ROPE[2],ROPE[2]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol)
     text( mean(ROPE), ROPEtextHt,
           bquote( .(round(100*pcInROPE))*"% in ROPE" ),
           adj=c(.5,0), cex=1, col=ropeCol )
  }
  # Display the HDI.
  HDI <- hdi( paramSampleVec, credMass )
  lines(HDI, c(0,0), lwd=4, lend='butt')
  text( mean(HDI), 0, bquote(.(100*credMass) * "% HDI" ),
        adj=c(.5,-1.7), cex=cex )
  text( HDI[1], 0, bquote(.(signif(HDI[1],3))),
        adj=c(HDItextPlace,-0.5), cex=cex )
  text( HDI[2], 0, bquote(.(signif(HDI[2],3))),
        adj=c(1.0-HDItextPlace,-0.5), cex=cex )

  return(invisible(NULL))
}
