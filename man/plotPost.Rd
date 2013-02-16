\name{plotPost}
\alias{plotPost}
\title{
Graphic display of a posterior probability distribution
}
\description{
Plot the posterior probability distribution for a single parameter from a vector of samples, typically from an MCMC process, with appropriate summary statistics.
}
\usage{
plotPost(paramSampleVec, credMass = 0.95, compVal = NULL, ROPE = NULL,
  HDItextPlace = 0.7, yaxt = "n", ylab = "", xlab = "Parameter", main = "", 
  cex.lab = 1.5, cex = 1.4, xlim = range(compVal, paramSampleVec),
  col = "skyblue", border = "white", showMode = FALSE, showCurve = FALSE,
  breaks = NULL, ...)
}
\arguments{
  \item{paramSampleVec}{
A vector of samples drawn from the target distribution.
}
  \item{credMass}{
the probability mass to include in credible intervals.
}
  \item{compVal}{
a value for comparison with those plotted.
}
  \item{ROPE}{
a two element vector, such as c(-1,1), specifying the limits of the Region Of Practical Equivalence. 
}
  \item{HDItextPlace}{
a value in [0,1] that controls the horizontal position of the labels at the ends of the HDI bar.
}
  \item{yaxt}{
a character which specifies the y axis type; the default, "n", suppresses plotting.
}
  \item{ylab}{
text to use as the label of the y axis.
}
  \item{xlab}{
text to use as the label of the x axis.
}
  \item{cex.lab}{
the magnification to be used for x and y labels relative to the current setting of cex
}
  \item{cex}{
a numerical value giving the amount by which plotting text and symbols should be magnified relative to the default
}
  \item{xlim}{
a vector of length 2 giving the limits for the x axis.
}
  \item{main}{
text to use as the main title of the plot
}
  \item{col}{
colour for the fill of the histogram boxes or for the curve if showCurve = TRUE.
}
  \item{border}{
colour for the borders of the histogram boxes; ignored if showCurve = TRUE.
}
  \item{showMode}{
logical: if TRUE, the mode is displayed instead of the mean.
}
  \item{showCurve}{
logical: if TRUE, the posterior density will be represented by a kernel density function instead of a histogram.
}
  \item{breaks}{
controls the histogram break points or the number of bars; see \code{\link{hist}}.
}
  \item{\dots}{
other graphical parameters
}
}
\details{
The data are plotted either as a histogram (above) or, if \code{showCurve = TRUE}, as a fitted kernel density curve (below). Either the mean or the mode of the distribution is displayed, depending on the parameter \code{showMode.} The Highest Density Interval (HDI) is shown as a horizontal bar, with labels for the ends of the interval.

\figure{plotPost1.jpg} \cr
\cr
\figure{plotPost2.jpg}

If values for a ROPE are supplied, these are shown as dark red vertical dashed lines, together with the percentage of probability mass within the ROPE. If a comparison value (\code{compVal}) is supplied, this is shown as a vertical green dotted line, together with the probability mass below and above this value.
}
\value{
Returns NULL invisibly. Used for its plotting side-effect.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
For details of the HDI calculation, see \code{\link{hdi}}.
}
\examples{
# Generate some data
tst <- rnorm(1e5, 3, 1)
plotPost(tst)
plotPost(tst, col='wheat', border='magenta')
plotPost(tst, credMass=0.8, ROPE=c(-1,1), xlab="Response variable")
plotPost(tst, showMode=TRUE, showCurve=TRUE, compVal=5.5)

}
