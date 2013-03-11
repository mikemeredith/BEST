\name{print.BEST}
\alias{print.BEST}
\title{
Printing a BEST object
}
\description{
Print method for objects of class BEST, such as produced by \code{\link{BESTmcmc}}.
}
\usage{
\method{print}{BEST}(x, digits = 4, ...)
}
\arguments{
  \item{x}{
an object of class \code{BEST}, as produced by \code{BESTmcmc}.
}
  \item{digits}{
the number of digits to print.
}
  \item{\dots}{
further arguments for the print function.
}
}
\details{
The print method displays summary statistics for the parameters and two MCMC diagnostic measures:

Rhat is the 'potential scale reduction factor', which is 1 on convergence; if any parameter has a value > 1.05, rerun with increased burnInSteps. See \code{\link[coda]{gelman.diag}}.

n.eff is the sample size adjusted for autocorrelation; for stable estimates of credible intervals this should be > 10,000. Rerun with increased numSavedSteps or increased thinSteps. See \code{\link[coda]{effectiveSize}}.
}
\value{
print.BEST returns x invisibly.
}
\author{
Mike Meredith
}

\seealso{
\code{\link{BESTmcmc}}.
}
\examples{
## See examples in BEST-package help.
}