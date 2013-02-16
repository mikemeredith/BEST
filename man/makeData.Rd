\name{makeData}
\alias{makeData}
\title{
Population parameter specification for a power analysis
}
\description{
The function allows the analyst to prepare an idealized data set which exactly matches selected point values, and incorporates uncertainty in these values in terms of sample size.
}
\usage{
makeData(mu1, sd1, mu2 = NULL, sd2 = NULL,
  nPerGrp, pcntOut = 0, sdOutMult = 2, rnd.seed = NULL)
}
\arguments{
  \item{mu1}{
the mean for the first (or only) population.
}
  \item{sd1}{
the standard deviation for the main part of the first population, excluding outliers.
}
  \item{mu2}{
the mean for the second population; NULL if only one population is involved.
}
  \item{sd2}{
the standard deviation for the main part of the second population; NULL if only one population is involved.
}
  \item{nPerGrp}{
sample size per group; large sample size reflects a high degree of precision in the values for the means and standard deviations.
}
  \item{pcntOut}{
the percentage of outliers in each population.
}
  \item{sdOutMult}{
the standard deviation of the outliers as a multiple of the standard deviation of the main part of the population.
}
  \item{rnd.seed}{
a seed for the random number generator, used to obtain reproducible samples if required.
}
}
\details{
The arguments to this function provide a framework to specify the hypothesised values of the parameters of the populations under study, while the sample size is chosen to reflect the confidence in the values specified.

The function produces idealized samples, ie. samples which exactly match the specified means and standard deviations, and displays the results as a plot.

\figure{makeData.jpg}

\emph{Histograms: actual sample values; red dashed line: distribution of the outliers; blue dashed line: distribution of the non-outliers; black line: combined distribution.}

These idealised samples are passed to \code{\link{BESTmcmc}}, which generates a series of sets of credible values for the parameters, including the normality parameter, taking account of correlations among them.

The sets of credible parameter values which constitute the \code{BESTmcmc} output are used by \code{\link{BESTpower}} to simulate new data sets which might arise during a subsequent experiment.
}
\value{
A list with two components:
 \item{y1}{A vector of simulated values for the first (or only) group.}
 \item{y2}{A vector of simulated values for the second group or NULL.}

}
\references{
Kruschke, J. K. 2012. Bayesian estimation supersedes the \emph{t} test. \emph{Journal of Experimental Psychology: General} First Posting. doi: 10.1037/a0029146
}
\author{
John Kruschke
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## See examples for BESTpower.
}
