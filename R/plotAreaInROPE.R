plotAreaInROPE <- 
function(paramSampleVec, credMass = 0.95, compVal = 0, maxROPEradius,
  plot = TRUE,...) {
  # Plots the probability mass included in the ROPE as a function of
  #   the half-width of the ROPE.

  # Sanity checks:
  if(missing(maxROPEradius))
    stop("maxROPEradius is missing with no default.")
  if(!isTRUE(maxROPEradius > 0))
    stop("maxROPEradius must be > 0.")
  if(!isTRUE(is.finite(compVal)))
    stop("A finite value for compVal is needed.")

  ropeRadVec = seq( 0 , maxROPEradius , length=201 ) # arbitrary comb
  areaInRope = rep( NA , length(ropeRadVec) )
  for ( rIdx in 1:length(ropeRadVec) ) {
    areaInRope[rIdx] <- mean( paramSampleVec > (compVal-ropeRadVec[rIdx])
                            & paramSampleVec < (compVal+ropeRadVec[rIdx]) )
  }

  if(plot) {
    plot( ropeRadVec , areaInRope ,
          xlab=bquote("Radius of ROPE around "*.(compVal)) ,
          ylab="Posterior in ROPE" ,
          type="l" , lwd=4 , col="darkred" , cex.lab=1.5 , ... )
    HDIlim = hdi( paramSampleVec , credMass=credMass )
    farHDIlim = HDIlim[which.max(abs(HDIlim-compVal))]
    ropeRadHDI = abs(compVal-farHDIlim)
    areaInFarHDIlim <- mean( paramSampleVec > (compVal-ropeRadHDI)
                           & paramSampleVec < (compVal+ropeRadHDI) )
    lines( c(ropeRadHDI, ropeRadHDI) , c(-0.5, areaInFarHDIlim) ,
           lty="dashed" , col="darkred" )
    text( ropeRadHDI , 0 ,
          bquote( atop( .(100*credMass)*"% HDI limit" ,
                       "farthest from "*.(compVal) ) ) , adj=c(0.5,0) )
    lines( c(-0.5, ropeRadHDI) ,c(areaInFarHDIlim, areaInFarHDIlim) ,
           lty="dashed" , col="darkred" )
    text( 0 , areaInFarHDIlim , bquote(.(signif(areaInFarHDIlim, 3))) ,
          adj=c(0, 1.1) )
  }

  invisible( list( x=ropeRadVec , y=areaInRope ) )
}
