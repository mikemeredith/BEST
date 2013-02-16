plot.BEST <-
function(x, which=c("mean", "sd", "effect", "nu"), credMass=0.95,
                    ROPE=NULL, compVal=0, showCurve=FALSE, ...) {
  # This function plots the posterior distribution for one selected item. 
  # Description of arguments:
  # x is mcmc.list object of the type returned by function BESTmcmc.
  # which indicates which item should be displayed; possible values are "mean", "sd",
  #  "effect" or "nu".
  # ROPE is a two element vector, such as c(-1,1), specifying the limit
  #   of the ROPE.
  # compVal is a scalar specifying the value for comparison.
  # showCurve if TRUE the posterior should be displayed as a fitted density curve
  #   instead of a histogram (default).

  # TODO additional sanity checks.

  xmat <- as.matrix(x)
  oneGrp <- ncol(xmat) == 3
  whichID <- match.arg(which)

  toPlot <- switch(whichID,
    mean = if(oneGrp) xmat[,"mu"] else xmat[,"mu[1]"] - xmat[,"mu[2]"],
    sd = if(oneGrp) xmat[,"sigma"] else xmat[,"sigma[1]"] - xmat[,"sigma[2]"],
    effect = if(oneGrp) (xmat[,"mu"] - compVal) / xmat[,"sigma"] else
       (xmat[,"mu[1]"] - xmat[,"mu[2]"]) /
          sqrt( ( xmat[,"sigma[1]"]^2 + xmat[,"sigma[2]"]^2 ) / 2 ),
    nu = log10(xmat[,"nu"]) )

  main <- switch(whichID,
    mean = if(oneGrp) "Mean" else "Difference of Means",
    sd = if(oneGrp) "Standard Deviation" else "Difference of Std. Dev.s",
    effect = "Effect Size",
    nu = "Normality")

  xlab <- switch(whichID,
    mean = if(oneGrp) bquote(mu) else bquote(mu[1] - mu[2]),
    sd = if(oneGrp) bquote(sigma) else bquote(sigma[1] - sigma[2]),
    effect = if(oneGrp) bquote( (mu-.(compVal)) / sigma ) else
      bquote( (mu[1]-mu[2]) / sqrt((sigma[1]^2 +sigma[2]^2 )/2 ) ),
    nu = bquote("log10("*nu*")"))

  if(whichID=="nu" && !is.null(compVal) && compVal == 0)
    compVal <- NULL
  if(whichID=="sd" && oneGrp && !is.null(compVal) && compVal == 0)
    compVal <- NULL

  # Plot posterior distribution of selected item:
  plotPost(toPlot, col="skyblue", credMass=credMass, showCurve=showCurve,
                  xlab=xlab , cex.lab = 1.75 , showMode=whichID != "mean",
                  compVal=compVal, main=main ) 

  return(invisible(NULL))
}
