\name{plotPostPred}
\alias{plotPostPred}
\title{
Plots for Posterior Predictive checks.
}
\description{
Plots a number (default 30) of credible t-distributions based on posterior values of the mean, standard deviation, and normality for each group, together with histograms of the data.
}
\usage{
plotPostPred(BESTobj, nCurvesToPlot = 30, mainColor="skyblue", dataColor="red")
}
\arguments{
  \item{BESTobj}{
an object of class \code{BEST}, as produced by the function \code{\link{BESTmcmc}}.
}
  \item{nCurvesToPlot}{
the number of posterior predictive curves to plot.
}
  \item{mainColor}{
an optional color name such as \code{"skyblue"} or a RGB specification such as \code{"#87CEEB"} that controls the color of the posterior prediction lines.
}
  \item{dataColor}{
an optional color name such as \code{"red"} or a RGB specification such as \code{"#FF0000"} that controls the color of the data histogram.
}
}

\value{
Nothing, used for its side effect.
}
\references{
Kruschke, J. K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146
}
\author{
John Kruschke, modified by Mike Meredith.
}

\examples{
## See examples in BEST-package.
}

\keyword{hplot}
