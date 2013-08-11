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
Version: \tab 0.1.4\cr
Date: \tab 2013-08-11\cr
License: \tab GPL-3\cr
}

The core function, \code{\link{BESTmcmc}}, generates posterior distributions to compare the means of two groups, or to compare the mean of one group with a standard, taking into account the standard deviation(s). It is thus similar to a \emph{t} test. However, our Bayesian approach results in probability statements about the values of interest, rather than \emph{p}-values and significance levels.

In addition, the procedure accounts for departures from normality by using a \emph{t}-distribution to model the variable of interest and estimating a measure of normality.

Functions to summarise and to visualise the output are provided.

The function \code{\link{BESTpower}} allows simulation-based estimates of power, either retrospective power directly with \code{BESTmcmc} output or prospective power analysis with \code{\link{makeData}}.

}
\author{
Original code by John K. Kruschke, \url{http://www.indiana.edu/~kruschke/BEST/}

Maintainer: Mike Meredith <mmeredith at wcs.org>
}

\references{
Kruschke, J. K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

Kruschke, J. K. 2011. \emph{Doing Bayesian data analysis: a tutorial with R and BUGS.} Elsevier, Amsterdam, especially Chapter 18.
}

\keyword{ package }
\keyword{ htest }

\examples{

\donttest{
## Comparison of two groups:
## =========================
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)

# Run an analysis, takes up to 1 min.
BESTout <- BESTmcmc(y1, y2)

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

head(BESTout$mu1)
muDiff <- BESTout$mu1 - BESTout$mu2
mean(muDiff > 1.5)
mean(BESTout$sigma1 - BESTout$sigma2)
hist(BESTout$nu)

# Retrospective power analysis
# ----------------------------
# This takes time, so we do 2 simulations here; a real analysis needs several hundred

powerRet <- BESTpower(BESTout, N1=length(y1), N2=length(y2), 
            ROPEm=c(-0.1,0.1), maxHDIWm=2.0, nRep=2) 
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

head(BESTout1$mu)
mean(BESTout1$sigma)
} % end of donttest
}
