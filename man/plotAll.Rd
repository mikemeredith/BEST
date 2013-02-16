\name{plotAll}
\alias{plotAll}
\title{
A series of plots displaying the results of 'BEST' analysis.
}
\description{
Displays a series of plots showing the posterior probability distributions of the parameters of interest.
}
\usage{
plotAll(BESTobj, y1=NULL, y2=NULL, credMass = 0.95,
  ROPEm = NULL, ROPEsd = NULL, ROPEeff = NULL,
  compValm = 0, compValsd = NULL, compValeff = 0,
  showCurve = FALSE, ...)
}
\arguments{
  \item{BESTobj}{
an object of class 'BEST', as produced by the function \code{\link{BESTmcmc}}.
}
  \item{y1}{
optional, the set of observations for the first (or only) group.
}
  \item{y2}{
optional, the set of observations for the second group, ignored if there is only one group.
}
  \item{credMass}{
the probability mass to include in credible intervals.
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
  \item{showCurve}{
logical: if TRUE, the posterior density will be represented by a kernel density function instead of a histogram.
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
  \item{...}{
other graphical parameters.
}
}
\details{
The display has a series of panels displaying the posterior distributions of each of the parameters (and differences between groups) together with summary statistics; see \code{\link{plotPost}} for details.
Also a chart showing approx. 30 plots of posterior predictive distributions, together with histograms of y1 and y2 if these are provided.
}
\value{
Returns NULL invisibly. Used for the side effect.
}
\references{
Kruschke, J. K. 2012. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} First Posting. doi: 10.1037/a0029146
}
\author{
Code by John Kruschke, modified by Mike Meredith.
}

\seealso{
\code{\link[=plot.BEST]{plot}} for plots of individual parameters. \code{\link[=summary.BEST]{summary}} for values of the corresponding summary statistics and \code{\link[=pairs.BEST]{pairs}} for a scatterplot matrix plot and correlation coefficients.
}
\examples{
# See examples in BEST-package.
}

