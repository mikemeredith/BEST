\name{plotAreaInROPE}
\alias{plotAreaInROPE}

\title{
Area of the posterior density in the ROPE as a function of its width.
}
\description{
Calculates and (optionally) plots the posterior probability mass included in the Region of Practical Equivalence (ROPE: see \code{\link{plot.BEST}}) as a function of the width of the ROPE.
}

\usage{
plotAreaInROPE(paramSampleVec, credMass = 0.95, compVal = 0, maxROPEradius,
  n = 201, plot = TRUE, ...)
}

\arguments{
  \item{paramSampleVec}{
A vector of samples drawn from the target distribution; see Examples.
}
  \item{credMass}{
The probability mass to include in credible intervals.
}
  \item{compVal}{
a value for comparison with those plotted.
}
  \item{maxROPEradius}{
The maximum value of the ROPE radius (ie. half-width) to include in the plot.
}
  \item{n}{
The number of equally spaced points at which the area in the ROPE is to be estimated.
}
  \item{plot}{
If FALSE, the plot will be suppressed but the values will be returned.
}
  \item{\dots}{
Other graphical parameters.
}
}

\details{
Defining a Region of Practical Equivalence (ROPE) allows decisions on whether a parameter is, for practical purposes, equivalent to a hypothetical null value, given a posterior probability density for the parameter. The null value may be considered credible if (A) 95\% (say) of the probability mass lies within the ROPE, or (B) the 95\% highest density interval (95\% HDI) lies entirely within the ROPE.

How wide should the ROPE be? Different people at different times will have different ideas on the range of values equivalent to the null. The function \code{plotAreaInROPE} plots the probability mass lying within the ROPE for a range of widths (or rather radii or half-widths). It also shows the radius at which the HDI falls entirely within the ROPE.
}
\value{
Returns invisibly a list with elements:
  \item{x }{A vector of ROPE radii from 0 to \code{maxROPEradius}.}
  \item{y }{The corresponding proportion of the posterior density included in the ROPE.}
}
\references{
\url{http://doingbayesiandataanalysis.blogspot.com/2013/08/how-much-of-bayesian-posterior.html}
}
\author{
John K. Kruschke, with minor modifications by Mike Meredith.
}

\examples{
# Generate a fake MCMC posterior for effect size and plot it:
mcmcChain <- rnorm(50000,0.03,0.025)
plotPost(mcmcChain, compVal=0, ROPE=c(-0.1, 0.1))

# How does the mass within the ROPE vary with ROPE radius?
plotAreaInROPE(mcmcChain, credMass = 0.95, compVal = 0,
  maxROPEradius = 0.15)

\donttest{
# Generate real MCMC chains, takes up to 1 min:
y1 <- c(4.77, 4.33, 3.59, 3.33, 2.66, 3.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
BESTout <- BESTmcmc(y1, y2)
plot(BESTout)

meanDiff <- BESTout$mu1 - BESTout$mu2
plotAreaInROPE(meanDiff, credMass = 0.95, compVal = 0,
  maxROPEradius = 3)
}
}
\keyword{hplot}
