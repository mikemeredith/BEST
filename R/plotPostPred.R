# Changed 11 Aug 2017 to show data as a typical histogram

plotPostPred <-
function(BESTobj, nCurvesToPlot = 30,
         mainColor="skyblue", dataColor="red") {
  # This function plots the posterior predictive distribution and the data.
  # Description of arguments:
  # BESTobj is mcmc.list object of the type returned by function BESTmcmc.

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


  # mcmcChain <- as.matrix(BESTobj)
  data <- attr(BESTobj, "data")

  # Set up window and layout:
  oldpar <- par(mar=c(3.5,3.5,2.5,0.5), mgp=c(2.25,0.7,0), "mfrow")
    on.exit(par(oldpar))
  if(!oneGrp)
    par(mfrow=2:1)

  # Select thinned steps in chain for plotting of posterior predictive curves:
  stepIdxVec <- seq(1, NROW( BESTobj ), length= nCurvesToPlot)
  toPlot <- BESTobj[stepIdxVec, ]

  plotDataPPC(toPlot=toPlot, oneGrp=oneGrp, data=data, lineColor = mainColor , dataColor = dataColor)

}
