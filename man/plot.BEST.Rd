\name{plot.BEST}
\alias{plot.BEST}
\title{
A plot method for objects of class 'BEST'
}
\description{
Displays a plot showing the posterior probability distribution of one of the parameters of interest, the difference in means by default.
}
\usage{
\method{plot}{BEST}(x, which = c("mean", "sd", "effect", "nu"), credMass = 0.95,
    ROPE = NULL, compVal = 0, showCurve = FALSE, ...)
}
\arguments{
  \item{x}{
an object of class \code{BEST}, as produced by the function \code{\link{BESTmcmc}}.
}
  \item{which}{
one of "mean", "sd", "effect" or "nu" or an abbreviation of one of these; indicates which estimate to plot. For a comparison of two groups, "mean" and "sd" display the difference in means or standard deviation.
}
   \item{credMass}{
the probability mass to include in credible intervals.
}
  \item{ROPE}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the estimate. 
}
  \item{showCurve}{
logical: if TRUE, the posterior density will be represented by a kernel density function instead of a histogram.
}
  \item{compVal}{
a value for comparison with the (difference of) parameters.
}
  \item{...}{
other graphical parameters.
}
}
\details{
The posterior distribution is shown as a histogram or density curve (if \code{showCurve = TRUE}), together with the Highest Density Interval. A ROPE and comparison value are also shown if appropriate.
}
\value{
Returns NULL invisibly. Used for the side effect.
}
\references{
Kruschke, J. K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146
}
\author{
Mike Meredith, adapted from code by John Kruschke.
}

\seealso{
\code{\link{plotAll}} for a single plot showing all the parameters, \code{\link[=summary.BEST]{summary}} for values of the corresponding summary statistics and \code{\link[=pairs.BEST]{pairs}} for a scatterplot matrix plot and correlation coefficients.
}
\examples{
# See examples in BEST-package.
}

