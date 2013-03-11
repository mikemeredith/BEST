\name{BEST-package}
\alias{BEST-package}
\alias{BEST}
\docType{package}
\title{
Bayesian Estimation Supersedes the t Test
}
\description{
An alternative to \emph{t} tests, producing posterior estimates for groups means and standard deviations and their differences and effect sizes. Bayesian estimation provides a much richer picture of the data, and can be summarised as point estimates and credible intervals.
}
\details{
\tabular{ll}{
Package: \tab BEST\cr
Type: \tab Package\cr
Version: \tab 0.02\cr
Date: \tab 2013-02-12\cr
License: \tab GPL-3\cr
}

The package provides a series of functions to compare the means of two groups, or to compare the mean of one group with a standard, taking into account the standard deviation(s). It is thus similar to a \emph{t} test. However, our Bayesian approach results in probability statements about the values of interest, rather than \emph{p}-values and significance levels.

In addition, the procedure accounts for departures from normality by using a \emph{t}-distribution to model the variable of interest and estimating a measure of normality.

The function \code{\link{BESTmcmc}} generates vectors of random draws from the posterior distributions of the group means (\eqn{\mu}) and standard deviations (\eqn{\sigma}), as well as the measure of normality (\eqn{\nu}). The procedure uses a Bayesian MCMC process implemented in JAGS (Plummer 2003) with vague priors for all parameters, as described by Kruschke (2012). The model is shown in the diagram below.

\figure{BESTmodel.jpg}

Derived parameters, including the differences in means or standard deviations, and effect sizes can be obtained from the results of the \code{BESTmcmc} run.

The output from \code{BESTmcmc} has class \code{BEST}, which has print, plot and summary methods. These permit the extraction and display of credible intervals and proportions of the posterior mass above or below values of interest.

Class \code{BEST} inherits from \code{\link[coda]{mcmc.list}}, so all the diagnostic functions in the \code{coda} package are available.
}
\author{
Original code by John K. Kruschke, johnkruschke@gmail.com

\url{http://www.indiana.edu/~kruschke/BEST/}

Maintainer: Mike Meredith <mmeredith@wcs.org>
}

\references{
Kruschke, J. K. 2012. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} First Posting. doi: 10.1037/a0029146

Kruschke, J. K. 2011. \emph{Doing Bayesian data analysis: a tutorial with R and BUGS.} Elsevier, Amsterdam

Plummer, Martyn (2003). JAGS: A Program for Analysis of Bayesian Graphical Models Using Gibbs Sampling, \emph{Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003)}, March 20-22, Vienna, Austria. ISSN 1609-395X
}

\keyword{ package }
\keyword{ htest }

\examples{
## Comparison of two groups:
## =========================
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)

# Run an analysis, takes up to 1 min.
BESTout <- BESTmcmc(y1, y2)

# coda diagnostics
gelman.diag(BESTout)
effectiveSize(BESTout)

# Look at the result:
BESTout
summary(BESTout)
plot(BESTout)
plot(BESTout, "sd")
plotPostPred(BESTout)
plotAll(BESTout, credMass=0.8, ROPEm=c(-0.1,0.1), 
          ROPEeff=c(-0.2,0.2), compValm=0.5) 
summary(BESTout, credMass=0.8, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
          ROPEeff=c(-0.2,0.2)) 
pairs(BESTout)

attachBEST(BESTout)
head(mu1)
muDiff <- mu1 - mu2
mean(muDiff > 1.5)
mean(sigma1 - sigma2)
hist(nu)

# Retrospective power analysis
# ----------------------------
# This takes time, so we do 2 simulations here; a real analysis needs several hundred

powerRet <- BESTpower(BESTout, N1=length(y1), N2=length(y2), 
            ROPEm=c(-0.1,0.1), maxHDIWm=2.0, 
            nRep=2, mcmcLength=10000, saveName = NULL) 
powerRet
# We only set criteria for the mean, so results for sd and effect size are all NA.

## Analysis with a single group:
## =============================
y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)

# Run an analysis, takes up to 40 secs.
BESTout1 <- BESTmcmc(y0)
BESTout1
summary(BESTout1)
plot(BESTout1)

attachBEST(BESTout1)
head(mu)
mean(sigma)
detachBEST()
}
