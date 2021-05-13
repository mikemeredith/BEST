\name{BESTmcmc}
\alias{BESTmcmc}
\title{
Generate MCMC samples for posterior distributions
}
\description{
This function is the core of the BEST package. It calls JAGS and passes a description of the model, priors, and data, then retrieves and returns the MCMC samples for the parameters.
}
\usage{

BESTmcmc(y1, y2 = NULL, priors = NULL, doPriorsOnly = FALSE,
  numSavedSteps = 1e+05, thinSteps = 1, burnInSteps = 1000,
	verbose=TRUE, rnd.seed=NULL, parallel=NULL)
}
\arguments{
  \item{y1}{
a numeric vector of data values.
}
  \item{y2}{
a vector of values for a second group, or NULL if there is only one group of observations.
}
 \item{priors}{
an optional list of values controlling the priors, see Details.
}
 \item{doPriorsOnly}{
if TRUE, \code{BESTmcmc} returns MCMC chains representing the prior distributions, \emph{not} the posterior distributions for your data set.
}
  \item{numSavedSteps}{
the number of MCMC observations to be returned.
}
  \item{thinSteps}{
thinning rate. If set to n > 1, n steps of the MCMC chain are calculated for each one returned. This is useful if autocorrelation is high and you need to run long chains.
}
  \item{burnInSteps}{
number of steps to discard as burn-in at the beginning of the chain.
}
   \item{verbose}{
if FALSE, output to the R Console is suppressed. If chains are run in parallel, the output from JAGS is not displayed in the Console, even if \code{verbose = TRUE}.
}
  \item{rnd.seed}{
a positive integer (or NULL): the seed for the random number generator, used to obtain reproducible samples if required. Values generated in different versions of BEST or different versions of JAGS may differ, even with the same seed.
}
  \item{parallel}{
if NULL or TRUE and > 3 cores are available, the MCMC chains are run in parallel. (If TRUE and < 4 cores are available, a warning is given.)
}
}
\details{
The function uses a t-distribution to model each sample, and generates vectors of random draws from the posterior distribution of the \emph{center} (\eqn{\mu}) and \emph{spread} or \emph{scale} (\eqn{\sigma}) of the distribution, as well as a measure of \emph{normality} (\eqn{\nu}).  The procedure uses a Bayesian MCMC process implemented in JAGS (Plummer 2003).

\eqn{\mu} is the population mean, except when \eqn{\nu} = 1 (which is the Cauchy distribution) or lower, when the mean is undefined.

\eqn{\sigma} is a good approximation to the standard deviation (SD) for values of \eqn{\nu} > 20. More exactly the SD is \eqn{\sigma} * sqrt(\eqn{\nu}/(\eqn{\nu} - 2)). For a normal distribution (with \eqn{\nu = \infty}), SD = \eqn{\sigma} is exact. The SD is undefined when \eqn{\nu} = 2 or less.

If \code{priors = NULL}, broad priors as described by Kruschke (2013) are used. For \eqn{\mu}, \code{Normal(mean(y), 1000 * sd(y))}; for \eqn{\sigma}, \code{Uniform(sd(y)/1000, sd(y) * 1000)}; for \eqn{\nu}, \code{Exponential(1/29) + 1}, with the constraint that \code{nu >= 1}. Here \code{y = c(y1, y2)}. Note that \code{priors = NULL} is not equivalent to \code{priors = list()}.

Alternatively, \code{priors} can be a list with elements specifying the priors for one or more parameters:\cr
\eqn{\mu} : population centers have separate normal priors, with mean \code{muM} and standard deviation \code{muSD}; if not included in the list, default values of \code{muM = mean(y), muSD = sd(y)*5} are used;\cr
\eqn{\sigma} : population scales have separate gamma priors, with \emph{mode} \code{sigmaMode} and standard deviation \code{sigmaSD}; defaults are \code{sigmaMode = sd(y), sigmaSD = sd(y)*5};\cr
\eqn{\nu} : the normality parameter has a gamma prior with \emph{mean} \code{nuMean} and standard deviation \code{nuSD}; defaults are \code{nuMean = 30, nuSD = 30}; versions before 0.4.0 constrained \eqn{\nu} to be >1.

If there are 2 groups of observations, \code{muM, muSD, sigmaMode, sigmaSD} may be vectors of length 2 or scalar; if scalar, the same value is used for each population.

The model is shown in the diagram below.

\figure{BESTmodel.jpg}

Derived parameters, including the differences in means or standard deviations, and effect sizes can be obtained from the results of the \code{BESTmcmc} run.

The output from \code{BESTmcmc} has class \code{BEST}, which has print, plot and summary methods. These permit the extraction and display of credible intervals and proportions of the posterior mass above or below values of interest.

}
\value{
An object of class \code{BEST} inheriting from \code{data.frame}. If two samples are compared, the output has the following columns:

  \item{mu1, mu2}{simulated observations of center for each population}
  \item{sigma1, sigma2}{simulated observations of scale for each population}
  \item{nu}{simulated observations of normality parameter}

while for a single sample, the columns are \code{mu, sigma, nu}.

The output has the following attributes:

  \item{call}{the call to the function.}
  \item{Rhat}{the 'potential scale reduction factor'.}
  \item{n.eff}{sample size adjusted for autocorrelation.}
  \item{data}{a list with elements y1 and y2 containing the original data; y2 may be NULL.}
  \item{priors}{a list with the priors used, if the \code{priors} argument is not NULL.}
  \item{doPriorsOnly}{logical, the value of the \code{doPriorsOnly} argument.}

The package provides \code{print}, \code{plot} and \code{summary} methods for \code{BEST} objects.
}
\references{
Kruschke, J K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

For the informative priors, see Kruschke's blog post at \url{http://doingbayesiandataanalysis.blogspot.com/2015/04/informed-priors-for-bayesian-comparison.html}

For the constraint on \eqn{\nu}, see the blog post at \url{http://doingbayesiandataanalysis.blogspot.com/2015/12/prior-on-df-normality-parameter-in-t.html}

Plummer, Martyn (2003). JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling, \emph{Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)}, March 20-22, Vienna, Austria. ISSN 1609-395X
}
\author{
Original code by John K. Kruschke, modified by Mike Meredith.
}

\seealso{
\code{\link[=plot.BEST]{plot}}, \code{\link[=summary.BEST]{summary}}, \code{\link[=pairs.BEST]{pairs}} for relevant methods.
}
\examples{
## See examples in BEST-package help.
}
