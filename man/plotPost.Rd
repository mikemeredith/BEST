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
  HDItextPlace = 0.7, showMode = FALSE, showCurve = FALSE, 
  mainColor="skyblue", comparisonColor="darkgreen", ROPEColor = "darkred",
  ...)
}
\arguments{
  \item{paramSampleVec}{
A vector of samples drawn from the target distribution.
}
  \item{credMass}{
the probability mass to include in credible intervals, or NULL to suppress plotting of credible intervals.
}
  \item{compVal}{
a value for comparison with those plotted.
}
  \item{ROPE}{
a two element vector, such as \code{c(-1, 1)}, specifying the limits of the Region Of Practical Equivalence. 
}
  \item{HDItextPlace}{
a value in [0,1] that controls the horizontal position of the labels at the ends of the HDI bar.
}
  \item{showMode}{
logical: if TRUE, the mode is displayed instead of the mean.
}
  \item{showCurve}{
logical: if TRUE, the posterior density will be represented by a kernel density function instead of a histogram.
}
  \item{mainColor}{
an optional color name such as \code{"skyblue"} or a RGB specification such as \code{"#87CEEB"} that controls the color of the bar charts and posterior prediction lines.
} 
  \item{comparisonColor}{
a optional color name such as \code{"darkgreen"} or a RGB specification such as \code{"#013220"} that controls the color of the text that specifies the relation of the HDI to a comparison value.
}  
  \item{ROPEColor}{
a optional color name such as \code{"darkred"} or a RGB specification such as \code{"#8B0000"} that controls the color of the text that specifies the relation of the HDI to the ROPE.
}
  \item{\dots}{
graphical parameters and the \code{breaks} parameter for the histogram.
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
Returns an object of class \code{histogram} invisibly. Used for its plotting side-effect.
}
\author{
John Kruschke, modified by Mike Meredith
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

# For integers:
tst <- rpois(1e5, 12)
plotPost(tst)

# A severely bimodal distribution:
tst2 <- c(rnorm(1e5), rnorm(5e4, 7))
plotPost(tst2)                  # A valid 95% CrI, but not HDI
plotPost(tst2, showCurve=TRUE)  # Correct 95% HDI

}
\keyword{hplot}