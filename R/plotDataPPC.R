# Function rewritten 11 Aug 2017 to deal with the whole process of doing 1 or 2 plots
#  with same xlim and ylim parameters.

plotDataPPC <-
function(toPlot, oneGrp, data) {
  # Does the plots of posterior predictive curves for one OR TWO samples
  # Called by plotAll and plotPostPred; no sanity checks; not exported.
  # Calling function should arrange for multiple plots with par(mfrow) or layout.
  # Now DOES do title and sample size.
  # toPlot : a data frame with parameters to use for the t-curves; for one sample case,
  #   colnames should be c("mu1","nu","sigma1")
  # oneGroup : TRUE for 1 group, FALSE for 2
  # data : list with components y1 and y2, y2 should be NULL in the one-sample case.

  # Work out the x axis limits
  if(is.null(data$y1) && is.null(data$y2)) {
    xRange <- range(toPlot$mu1, toPlot$mu2)
  } else {
    # Get the breaks for the histograms, both must be same
    breaks <- hist(c(data$y1, data$y2), plot=FALSE)$breaks
    xRange <- range(breaks)
  }
  xLim <- c( xRange[1]-0.1*diff(xRange) ,
              xRange[2]+0.1*diff(xRange) )

  # Prepare the stuff to plot, so we can get the y axis limit
  npoints <- 100
  nlines <- nrow(toPlot)
  nplots <- if(oneGrp) 1 else 2
  xVec <- seq(xLim[1], xLim[2], length=npoints)
  PPDmat <- array(NA, c(npoints, nlines, nplots))
  for(i in 1:nlines)
    PPDmat[, i, 1] <- dt( (xVec-toPlot$mu1[i])/toPlot$sigma1[i], df=toPlot$nu[i] )/toPlot$sigma1[i]
  if(!oneGrp)
    for(i in 1:nlines)
      PPDmat[, i, 2] <- dt( (xVec-toPlot$mu2[i])/toPlot$sigma2[i], df=toPlot$nu[i] )/toPlot$sigma2[i]
  hist1 <- hist2 <- NULL
  if(!is.null(data$y1))
    hist1 <- hist(data$y1, breaks=breaks, plot=FALSE)
  if(!is.null(data$y2))
    hist2 <- hist(data$y2, breaks=breaks, plot=FALSE)
  # Now get y axis limit
  maxY <- max(PPDmat, hist1$density, hist2$density)

  # Do first plot
  plot(xVec[1], 0, xlim=range(xVec), ylim=c(0, maxY), cex.lab=1.75,
        type="n", xlab="y", ylab="p(y)", lwd=1)
  if(oneGrp) {
    title(main="Data w. Post. Pred.")
    if(!is.null(data$y1))
      text( max(xVec) , maxY , bquote(N ==.(length(data$y1))) , adj=c(1.1,1.1) )
  } else {
    title(main="Data Group 1 w. Post. Pred.")
    if(!is.null(data$y1))
      text( max(xVec) , maxY , bquote(N[1]==.(length(data$y1))) , adj=c(1.1,1.1) )
  }
  matlines(x=xVec, y=PPDmat[, , 1], lty=1, col="skyblue")
  if(!is.null(hist1)) {
    op <- par(lwd=2)
    plot(hist1, freq=FALSE, border='red', add=TRUE)
    segments(x0=xVec[1], y0=0, x1=xVec[npoints], col='red')
    par(op)
  }
  # Maybe do second plot
  if(!oneGrp) {
    plot(xVec[1], 0, xlim=range(xVec), ylim=c(0, maxY), cex.lab=1.75,
          type="n", xlab="y", ylab="p(y)")
    title(main="Data Group 2 w. Post. Pred.")
    if(!is.null(data$y2))
      text( max(xVec) , maxY , bquote(N[2]==.(length(data$y2))) , adj=c(1.1,1.1) )

    matlines(x=xVec, y=PPDmat[, , 2], lty=1, col="skyblue")
    if(!is.null(hist2)) {
      op <- par(lwd=2)
      plot(hist2, freq=FALSE, border='red', add=TRUE)
      segments(x0=xVec[1], y0=0, x1=xVec[npoints], col='red')
      par(op)
    }
  }
}
