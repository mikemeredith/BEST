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
 saveName = "BESTpower.Rdata", showFirstNrep = 0)
}
\arguments{
  \item{BESTobj}{
an object of class BEST produced by BESTmcmc.
}
  \item{N1}{
planned sample size for the first (or only) group of observations; may be a scalar if sample size is fixed, or a vector if sample size varies; values will be recycled if necessary.
}
  \item{N2}{
planned sample size for the second group of observations; ignored if BESTobj concerns only one group.
}
  \item{credMass}{
the probability mass to include in HDIs when checking criteria.  
}
  \item{ROPEm}{
a two element vector, such as c(-1,1), specifying the limit of the ROPE on the difference of means (for 2 groups) or the mean (for 1 group). 
}
  \item{ROPEsd}{
a two element vector, such as c(-1,1), specifying the limit of the ROPE on the (difference of) standard deviations.
}
  \item{ROPEeff}{
a two element vector, such as c(-1,1), specifying the limit of the ROPE on the effect size.
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
by default, the results are saved to a file after each iteration; \code{saveName} specifies the file name to use. Set to NULL to disable saving.
}
  \item{showFirstNrep}{
the number of results to display as plots at the beginning of the simulation run. (Note for Mac users: this uses X11 graphics which is no longer included in Mac OS X 10.8. If you try to run it, you will be prompted to install XQuartz, which should then work.) 
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
At least 1000 simulations are needed to get good estimates of power and these can take a long time. If the run is interrupted, the results so far can be recovered from the file specified in \code{saveName}.

The chains in \code{BESTobj} must have at least nRep values. To allow for some degree of autocorrelation among values, it would be prudent to make these chains at least 10 * nRep in length.
}

\section{To Do}{
1. It would be nice to program this to run on parallel cores.

2. Check that the plots work on platforms other than Windows.
}

\seealso{
\code{\link{makeData}} for details of preparing a \code{BESTobj} for a prospective power analysis.
}
\examples{

## For retrospective power analysis, see the example in BEST-package.

# 1. Generate idealised data set:
proData <- makeData(mu1=108, sd1=17, mu2=100, sd2=15, nPerGrp=20, 
                         pcntOut=10, sdOutMult=2.0, rnd.seed=NULL)

# 2. Generate credible parameter values from the idealised data:
proMCMC <- BESTmcmc(proData$y1, proData$y2, numSavedSteps=2000)  

# 3. Compute the prospective power for planned sample sizes:
# We'll  do just 5 simulations to show it works; should be several hundred.
N1plan <- N2plan <- 50
powerPro <- BESTpower(proMCMC, N1=N1plan, N2=N2plan,
               ROPEm=c(-1.5,1.5) , ROPEsd=c(-2,2) , ROPEeff=c(-0.5,0.5) , 
               maxHDIWm=15.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 ,
               nRep=5, saveName=NULL)
powerPro

}
