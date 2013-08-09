\name{BESTpower}
\alias{BESTpower}
\title{
Estimating statistical power 
}
\description{
Estimation of the probability of meeting the goals of a study given initial information or assumptions about the population parameters. For prospective power estimation, the sequence\cr \code{makeData -> BESTmcmc -> BESTpower} \cr is recommended: see \code{\link{makeData}}.
}
\usage{
BESTpower(BESTobj, N1, N2, credMass=0.95,
 ROPEm, ROPEsd, ROPEeff, 
 maxHDIWm, maxHDIWsd, maxHDIWeff,
 compValm = 0, nRep = 200, mcmcLength = 10000,
 saveName = NULL, showFirstNrep = 0, verbose = 2, rnd.seed=NULL)
}
\arguments{
  \item{BESTobj}{
an object of class \code{BEST} produced by \code{BESTmcmc}.
}
  \item{N1}{
planned sample size for the first (or only) group of observations; may be a scalar if sample size is fixed, or a vector if sample size varies; values will be recycled if necessary.
}
  \item{N2}{
planned sample size for the second group of observations; ignored if \code{BESTobj} concerns only one group.
}
  \item{credMass}{
the probability mass to include in HDIs when checking criteria.  
}
  \item{ROPEm}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the difference of means (for 2 groups) or the mean (for 1 group). 
}
  \item{ROPEsd}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the (difference of) standard deviations.
}
  \item{ROPEeff}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the effect size.
}
  \item{maxHDIWm}{
the maximum acceptable width for the HDI for the difference in means (for 2 groups) or for the mean (for a single group).
}
  \item{maxHDIWsd}{
the maximum acceptable width for the HDI for the (difference of) standard deviation.
}
  \item{maxHDIWeff}{
the maximum acceptable width for the HDI for the effect size.
}
  \item{compValm}{
for a single group, the value of the mean which represents no effect; used to calculate the effect size. Ignored for 2 groups.
}
  \item{nRep}{
number of simulations to carry out.
}
  \item{mcmcLength}{
length of the MCMC chains to use for each simulation.
}
  \item{saveName}{
if required, the results may saved to a file after each iteration and \code{saveName} specifies the file name (or path relative to the current working directory) to use. The \code{power} object can be loaded with \code{\link{load}}. Set to NULL (the default) to disable saving.
}
  \item{showFirstNrep}{
the number of results to display as plots at the beginning of the simulation run. (This uses dev.new(), which does not work in Rstudio. The plots will appear sequentially in the plot window and you will have to use the back arrow to review them.) 
}
  \item{verbose}{
controls output to the R Console: 0 suppresses all output; 1 gives just a progress bar; 2 gives maximum detail.
}
  \item{rnd.seed}{
a positive integer (or NULL): the seed for the random number generator, used to obtain reproducible samples if required.
}
}
\details{
For each of the parameters of interest - (difference in) mean, (difference in) standard deviation and effect size - we consider 4 criteria and the probability that each will be met:

1. The HDI of the posterior density of the parameter lies entirely outside the ROPE and is greater than the ROPE. 

2. The HDI of the posterior density of the parameter lies entirely outside the ROPE and is less than the ROPE. 

3. The HDI of the posterior density of the parameter lies entirely inside the ROPE. 

4. The width of the HDI is less than the specified \code{maxHDIWx}.

The mass inside the above HDIs depends on the \code{credMass} argument.

A uniform beta prior is used for each of these probabilities and combined with the results of the simulations to give a conjugate beta posterior distribution. The means and 95\% HDI credible intervals are returned. 
}
\value{
A matrix with a row for each criterion and columns for the mean and lower and upper limits of a 95\% credible interval for the posterior probability of meeting the criterion.

Note that this matrix always has 12 rows. Rows corresponding to criteria which are not specified will have NAs.
}
\references{
Kruschke, J. K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

Kruschke, J. K. 2011. \emph{Doing Bayesian data analysis: a tutorial with R and BUGS.} Elsevier, Amsterdam, Chapter 13.
}
\author{
Original code by John Kruschke, modified by Mike Meredith.
}
\note{
At least 1000 simulations are needed to get good estimates of power and these can take a long time. If the run is interrupted, the results so far can be recovered from the file specified in \code{saveName}.

The chains in \code{BESTobj} must have at least nRep values. To allow for some degree of autocorrelation among values, it would be prudent to make these chains at least 10 * nRep in length.
}

\seealso{
\code{\link{makeData}} for details of preparing a \code{BESTobj} for a prospective power analysis.
}
\examples{

## For retrospective power analysis, see the example in BEST-package.

# 1. Generate idealised data set:
proData <- makeData(mu1=108, sd1=17, mu2=100, sd2=15, nPerGrp=20, 
                         pcntOut=10, sdOutMult=2.0, rnd.seed=NULL)
\donttest{
# 2. Generate credible parameter values from the idealised data:
proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=2000)  

# 3. Compute the prospective power for planned sample sizes:
# We'll  do just 5 simulations to show it works; should be several hundred.
N1plan <- N2plan <- 50
powerPro <- BESTpower(proMCMC, N1=N1plan, N2=N2plan,
               ROPEm=c(-1.5,1.5), ROPEsd=c(-2,2), ROPEeff=c(-0.5,0.5), 
               maxHDIWm=15.0, maxHDIWsd=10.0, maxHDIWeff=1.0, nRep=5)
powerPro
} ## End donttest
}
