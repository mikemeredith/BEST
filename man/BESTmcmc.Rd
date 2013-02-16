\name{BESTmcmc}
\alias{BESTmcmc}
\title{
Generate MCMC samples for posterior distributions
}
\description{
This function is the core of the BEST package. It calls JAGS and passes a description of the model and the data, then retrieves and returns the MCMC samples for the parameters.
}
\usage{
BESTmcmc(y1, y2 = NULL, numSavedSteps = 1e+05, thinSteps = 1, burnInSteps = 1000)
}
\arguments{
  \item{y1}{
a (non-empty) numeric vector of data values.
}
  \item{y2}{
a vector of values for a second group, or NULL if there is only one group of observations.
}
  \item{numSavedSteps}{
the number of MCMC observations to be returned.
}
  \item{thinSteps}{
thinning rate. If set to n > 1, one value in n steps of the MCMC chain is returned. This is useful if autocorrelation is high.
}
  \item{burnInSteps}{
number of steps to discard as burn-in at the beginning of the chain.
}
}
\details{

}
\value{
An object of class 'BEST' inheriting from mcmc.list.
The package provides print, plot and summary methods for 'BEST' objects.
}
\references{
Kruschke, J K. 2012. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} First Posting. doi: 10.1037/a0029146
}
\author{
Original code by John K. Kruschke, mangled by Mike Meredith
}
\note{
This differs from the original BESTmcmc function in Kruschke's code:\cr
1. Both 1-group and 2-group comparisons are dealt with by the same function.\cr
2. The output inherits from an mcmc.list, so functions in package coda can be used.
}


\seealso{
\code{\link[=plot.BEST]{plot}}, \code{\link[=summary.BEST]{summary}}, \code{\link[=pairs.BEST]{pairs}} for relevant methods.
}
\examples{
## See examples in BEST-package help.
}
