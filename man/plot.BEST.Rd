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
    ROPE = NULL, compVal = 0, showCurve = FALSE, 
    mainColor="skyblue", dataColor="red", comparisonColor="darkgreen", ROPEColor = "darkred",
    ...)
}
\arguments{
  \item{x}{
an object of class \code{BEST}, as produced by the function \code{\link{BESTmcmc}}.
}
  \item{which}{
one of "mean", "sd", "effect" or "nu" or an abbreviation of one of these; indicates which estimate to plot. For a comparison of two groups, "mean" and "sd" display the difference in means or standard deviation.
}
   \item{credMass}{
the probability mass to include in credible intervals, or NULL to suppress plotting of the credible interval.
}
  \item{ROPE}{
a two element vector, such as \code{c(-1, 1)}, specifying the limit of the ROPE on the estimate; see Details.
}
  \item{showCurve}{
logical: if TRUE, the posterior density will be represented by a kernel density function instead of a histogram.
}
  \item{compVal}{
a value for comparison with the (difference of) parameters.
}
  \item{mainColor}{
an optional color name such as \code{"skyblue"} or a RGB specification such as \code{"#87CEEB"} that controls the color of the bar charts and posterior prediction lines.
}
  \item{comparisonColor}{
a optional color name such as \code{"darkgreen"} or a RGB specification such as \code{"#013220"} that controls the color of the text that specifies the relation of the HDI to a comparison value.
}  
  \item{ROPEColor}{
a optional color name such as \code{"darkred"} or a RGB specification such as \code{"#8B0000"} that controls the color of the text that specifies the relation of the HDI to the ROPE.
} 
  \item{...}{
other graphical parameters.
}
}
\details{
The posterior distribution is shown as a histogram or density curve (if \code{showCurve = TRUE}), together with the Highest Density Interval. A ROPE and comparison value are also shown if appropriate.

The probability that the difference in means (or effect size, etc.) is precisely zero is zero. More interesting is the probability that the difference may be too small to matter. We can define a region of practical equivalence (ROPE) around zero, and obtain the posterior probability that the true value lies therein.
}
\value{
Returns an object of class \code{histogram} invisibly. Used mainly for the side effect.
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

