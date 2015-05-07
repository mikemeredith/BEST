\name{BESTmcmc}
\alias{BESTmcmc}
\title{
Generate MCMC samples for posterior distributions
}
\description{
This function is the core of the BEST package. It calls JAGS and passes a description of the model, priors, and data, then retrieves and returns the MCMC samples for the parameters.
}
\usage{
  
BESTmcmc(y1, y2 = NULL, priors = NULL, 
  numSavedSteps = 1e+05, thinSteps = 1, burnInSteps = 1000,
	verbose=TRUE, rnd.seed=NULL, parallel=TRUE)
}
\arguments{
  \item{y1}{
a numeric vector of at least 3 data values.
}
  \item{y2}{
a vector of at least 3 values for a second group, or NULL if there is only one group of observations.
}
 \item{priors}{
an optional list of values controlling the priors, see Details.
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
   \item{verbose}{
if FALSE, output to the R Console is suppressed.
}
  \item{rnd.seed}{
a positive integer (or NULL): the seed for the random number generator, used to obtain reproducible samples if required. Values generated in different versions of BEST may differ.
}
  \item{parallel}{
if TRUE and > 3 cores are available, the MCMC chains are run in parallel.
}
}
\details{
The function generates vectors of random draws from the posterior distributions of the group means (\eqn{\mu}) and standard deviations (\eqn{\sigma}), as well as the measure of normality (\eqn{\nu}). The procedure uses a Bayesian MCMC process implemented in JAGS (Plummer 2003).

BESTmcmc allows the user to specify priors for all parameters by including elements in the \code{priors} list:\cr
\eqn{\mu} : population means have separate normal priors, with mean \code{muM} and standard deviation \code{muSD};\cr
\eqn{\sigma} : population standard deviations have separate gamma priors, with \emph{mode} \code{sigmaMode} and standard deviation \code{sigmaSD};\cr
\eqn{\nu} : the normality parameter must be >= 1, so a gamma prior is used for \code{nu - 1}, with \emph{mean} \code{nuMean} and standard deviation \code{nuSD}.

If there are 2 groups of observations, \code{muM, muSD, sigmaMode, sigmaSD} may be vectors of length 2 or scalar, in which case the same value is used for each population. The default priors correspond to \code{list(muM = mean(y), muSD = sd(y)*5, sigmaMode = sd(y), sigmaSD = sd(y)*5, nuMean = 29, nuSD = 29)}, where \code{y = c(y1, y2)}.

If \code{priors = NULL}, vague priors as described by Kruschke (2013) are used for all parameters.


The model is shown in the diagram below.

\figure{BESTmodel.jpg}

Derived parameters, including the differences in means or standard deviations, and effect sizes can be obtained from the results of the \code{BESTmcmc} run.

The output from \code{BESTmcmc} has class \code{BEST}, which has print, plot and summary methods. These permit the extraction and display of credible intervals and proportions of the posterior mass above or below values of interest.

}
\value{
An object of class \code{BEST} inheriting from \code{data.frame}. If two samples are compared, the output has the following columns:

  \item{mu1, mu2}{simulated observations of means for each population}
  \item{sigma1, sigma2}{simulated observations of standard deviations for each population}
  \item{nu}{simulated observations of normality}

while for a single sample, the columns are \code{mu, sigma, nu}.

The output has the following attributes:

  \item{Rhat}{the 'potential scale reduction factor'}
  \item{n.eff}{sample size adjusted for autocorrelation}
  \item{data}{a list with elements y1 and y2 containing the original data; y2 may be NULL}
  \item{priors}{a list with the priors used, if the \code{priors} argument is not NULL}

The package provides \code{print}, \code{plot} and \code{summary} methods for \code{BEST} objects.
}
\references{
Kruschke, J K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

Plummer, Martyn (2003). JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling, \emph{Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)}, March 20-22, Vienna, Austria. ISSN 1609-395X
}
\author{
Original code by John K. Kruschke, modified by Mike Meredith.
}
\note{
This differs from the original \code{BESTmcmc} function in Kruschke's code:\cr
1. Both 1-group and 2-group comparisons are dealt with by the same function.\cr
2. The output inherits from \code{data.frame} rather then being a matrix. 
}

\seealso{
\code{\link[=plot.BEST]{plot}}, \code{\link[=summary.BEST]{summary}}, \code{\link[=pairs.BEST]{pairs}} for relevant methods.
}
\examples{
## See examples in BEST-package help.
}
