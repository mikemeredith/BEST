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
The function generates vectors of random draws from the posterior distributions of the group means (\eqn{\mu}) and standard deviations (\eqn{\sigma}), as well as the measure of normality (\eqn{\nu}). The procedure uses a Bayesian MCMC process implemented in JAGS (Plummer 2003) with vague priors for all parameters, as described by Kruschke (2012). The model is shown in the diagram below.

\figure{BESTmodel.jpg}

Derived parameters, including the differences in means or standard deviations, and effect sizes can be obtained from the results of the \code{BESTmcmc} run.

The output from \code{BESTmcmc} has class \code{BEST}, which has print, plot and summary methods. These permit the extraction and display of credible intervals and proportions of the posterior mass above or below values of interest.

Class \code{BEST} inherits from \code{\link[coda]{mcmc.list}}, so all the diagnostic functions in the \code{coda} package are available.
}
\value{
An object of class \code{BEST} inheriting from \code{mcmc.list} with the following additional attributes:

  \item{Rhat}{the 'potential scale reduction factor'}
  \item{n.eff}{sample size adjusted for autocorrelation}
  \item{data}{a list with elements y1 and y2 containing the original data; y2 may be NULL}

The package provides \code{print}, \code{plot} and \code{summary} methods for \code{BEST} objects.
}
\references{
Kruschke, J K. 2012. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} First Posting. doi: 10.1037/a0029146

Plummer, Martyn (2003). JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling, \emph{Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)}, March 20-22, Vienna, Austria. ISSN 1609-395X
}
\author{
Original code by John K. Kruschke, modified by Mike Meredith
}
\note{
This differs from the original \code{BESTmcmc} function in Kruschke's code:\cr
1. Both 1-group and 2-group comparisons are dealt with by the same function.\cr
2. In the present version, the output inherits from \code{mcmc.list}. This may change in future versions.
}

\seealso{
\code{\link[=plot.BEST]{plot}}, \code{\link[=summary.BEST]{summary}}, \code{\link[=pairs.BEST]{pairs}} for relevant methods.
}
\examples{
## See examples in BEST-package help.
}
