\name{summary.BEST}
\alias{summary.BEST}

\title{
Extract summary statistics from an object of class BEST.
}
\description{
Provides summary statistics for each of the parameters (mean and standard deviation) of the group(s) of observations and their differences.
}
\usage{
\method{summary}{BEST}(object, credMass = 0.95,
  ROPEm = NULL, ROPEsd = NULL, ROPEeff = NULL,
  compValm = 0, compValsd = NULL, compValeff = 0, ...)

}
\arguments{
  \item{object}{
an object of class \code{BEST}, as produced by the function \code{\link{BESTmcmc}}.
}
  \item{credMass}{
the probability mass to include in credible intervals.
}
  \item{ROPEm}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the difference of means (for 2 groups) or the mean (for 1 group). See \code{\link{plot.BEST}} for an explanation of ROPE.
}
  \item{ROPEsd}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the (difference of) standard deviations.
}
  \item{ROPEeff}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the effect size.
}
  \item{compValm}{
a value for comparison with the (difference of) means.
}
  \item{compValsd}{
a value for comparison with the (difference of) standard deviations.
}
  \item{compValeff}{
a value for comparison with the effect size.
}
  \item{\dots}{
additional arguments for the summary or print function.
}

}
\value{
Returns a matrix with the parameters in rows and the following columns:
\item{mean, median, mode}{ the mean, median and mode of the MCMC samples for the corresponding parameter.}
\item{hdi\%, hdiLow, hdiHigh}{the percentage of posterior probability mass included in the highest density interval and the lower and upper limits.}
\item{compVal, \%>compVal}{the value for comparison and the percentage of the posterior probability mass above that value.}
\item{ROPElow, ROPEhigh, \%InROPE}{the lower and upper limits of the Region Of Practical Equivalence (ROPE) and the percentage of the posterior probability mass within the region.}

If the analysis concerns a comparison of two groups, the matrix will have rows for:
\item{mu1, mu2, muDiff}{the means of each group and the difference in means}
\item{sigma1, sigma2, sigmaDiff}{the standard deviations of each group and the difference in standard deviations}
\item{nu, log10nu}{the normality parameter and its log}
\item{effSz}{the effect size; \eqn{d[a]} from Macmillan & Creelman (1991).}

For a single group, the rows will be:
\item{mu}{the mean}
\item{sigma}{the standard deviation}
\item{nu, log10nu}{the normality parameter and its log}
\item{effSz}{the effect size.}

Many of the elements of the matrix will be NA. The print method for the summary attempts to print this nicely.
}
\references{
Kruschke, J K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146

Macmillan, N. A., & Creelman, C. D. (1991). \emph{Detection Theory: A User's Guide}. New York, Cambridge University Press
}
\author{
Mike Meredith, based on code by John K. Kruschke.
}

\seealso{
Use the \code{\link{plotAll}} function for a graphical display of these same values.
}
\examples{
## see "BEST-package"
}

\keyword{methods}