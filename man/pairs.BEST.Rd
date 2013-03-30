\name{pairs.BEST}
\alias{pairs.BEST}
\title{
Scatterplot matrix for a \code{BEST} object
}
\description{
Function to produce a scatterplot matrix of a \code{BEST} object produced by \code{\link{BESTmcmc}}, with correlation coefficients in the lower triangle.
}
\usage{
\method{pairs}{BEST}(x, nPtToPlot = 1000, col = "skyblue", ...)
}
\arguments{
  \item{x}{
an object of class \code{BEST}
}
  \item{nPtToPlot}{
number of points to plot
}
  \item{col}{
colour to use for the points plotted.
}
  \item{\dots}{
other graphical parameters passed to \code{plot.default}.
}
}

\value{
None; used for its side effect.
}

\author{
Original code by John Kruschke, adapted as a \code{pairs} method by Mike Meredith
}

\seealso{
\code{\link{pairs}} in package \code{graphics}.
}
\examples{
# See examples in BEST-package
}
