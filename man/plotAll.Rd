\name{plotAll}
\alias{plotAll}
\title{
A series of plots displaying the results of 'BEST' analysis.
}
\description{
Displays a series of plots showing the posterior probability distributions of the parameters of interest.
}
\usage{
plotAll(BESTobj, credMass = 0.95,
  ROPEm = NULL, ROPEsd = NULL, ROPEeff = NULL,
  compValm = 0, compValsd = NULL, compValeff = 0,
  showCurve = FALSE, 
  mainColor="skyblue", dataColor="red", comparisonColor="darkgreen", ROPEColor = "darkred",...)
}
\arguments{
  \item{BESTobj}{
an object of class \code{BEST}, as produced by the function \code{\link{BESTmcmc}}.
}
  \item{credMass}{
the probability mass to include in credible intervals, or NULL to suppress plotting of the credible interval.
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
  \item{mainColor}{
an optional color name such as \code{"skyblue"} or a RGB specification such as \code{"#87CEEB"} that controls the color of the bar charts and posterior prediction lines.
} 
  \item{dataColor}{
an optional color name such as \code{"red"} or a RGB specification such as \code{"#FF0000"} that controls the color of the data histogram.
}  
  \item{comparisonColor}{
an optional color name such as \code{"darkgreen"} or a RGB specification such as \code{"#013220"} that controls the color of the text that specifies the relation of the HDI to a comparison value.
}  
  \item{ROPEColor}{
an optional color name such as \code{"darkred"} or a RGB specification such as \code{"#8B0000"} that controls the color of the text that specifies the relation of the HDI to the ROPE.
}  
  \item{...}{
other graphical parameters (currently ignored).
}
}
\details{
The display has a series of panels displaying the posterior distributions of each of the parameters (and differences between groups) together with summary statistics; see \code{\link{plotPost}} for details.
Also a chart showing approx. 30 plots of posterior predictive distributions, together with histograms of the original data.
}
\value{
Returns NULL invisibly. Used for the side effect.
}
\references{
Kruschke, J. K. 2013. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} 142(2):573-603. doi: 10.1037/a0029146
}
\author{
Code by John Kruschke, modified by Mike Meredith.
}

\seealso{
\code{\link[=plot.BEST]{plot}} for plots of individual parameters,  \code{\link[=summary.BEST]{summary}} for values of the corresponding summary statistics and \code{\link[=pairs.BEST]{pairs}} for a scatterplot matrix plot and correlation coefficients.
}
\examples{
# See examples in BEST-package.
}

