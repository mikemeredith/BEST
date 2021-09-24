\name{postPriorOverlap}
\alias{postPriorOverlap}
\title{
Overlap between posterior and prior probability distributions.
}
\description{
Calculates and displays the overlap between a posterior distribution (as a vector of samples, typically from an MCMC process) and a prior distribution (as a vector of samples or as a function). Unidentifiable parameters will have high overlap: Gimenez et al (2009) suggest that overlap greater than 35\% indicates weak identifiability.
}
\usage{
postPriorOverlap(paramSampleVec, prior, ..., yaxt="n", ylab="",
           xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
           xlim=range(paramSampleVec), breaks=NULL,
           mainColor="skyblue", priorColor="yellow", overlapColor="green")
}
\arguments{
  \item{paramSampleVec}{
a vector of samples drawn from the target distribution.
}
  \item{prior}{
\emph{either} a vector of samples drawn from the prior distribution \emph{or} the name for the density function of the distribution; standard R functions for this have a \code{d-} prefix, eg. \code{dbeta}. Arguments required by the function must be specified by their (abbreviated) names in the \code{\dots} argument; see the examples.
}
  \item{...}{
named parameters to be passed to \code{prior} when it is a function.
}
  \item{yaxt}{
a character which specifies the y axis type; the default, "n", suppresses plotting.
}
  \item{ylab}{
text to use as the label of the y axis.
}
  \item{xlab}{
text to use as the label of the x axis.
}
  \item{cex.lab}{
the magnification to be used for x and y labels relative to the current setting of \code{cex}
}
  \item{cex}{
a numerical value giving the amount by which plotting text and symbols should be magnified relative to the default
}
  \item{xlim}{
a vector of length 2 giving the limits for the x axis.
}
  \item{main}{
text to use as the main title of the plot
}
  \item{breaks}{
controls the histogram break points or the number of bars; see \code{\link{hist}}.
}
  \item{mainColor}{
an optional color name such as \code{"skyblue"} or a RGB specification such as \code{"#87CEEB"} that controls the color of the bar chart representing the posterior.
}
  \item{priorColor}{
an optional color name such as \code{"yellow"} or a RGB specification such as \code{"#FFFF00"} that controls the color of prior, both if it is data and when it is a function.
}
  \item{overlapColor}{
an optional color name such as \code{"green"} or a RGB specification such as \code{"#00FF00"} that controls the color of the overlap area.
} 
}
\value{
Returns the overlap, the area lying under the lower of the two density curves.
}
\references{
Gimenez, Morgan and Brooks (2009) Weak identifiability in models for mark-recapture-recovery data. pp.1055-1068 in Thomson, Cooch and Conroy (eds) \emph{Modeling demographic processes in marked populations} Springer
}
\author{
Mike Meredith
}

\examples{
# Generate some data
tst <- rbeta(1e6, 5, 7)

# check overlap with a Beta(0.2, 0.2) prior:
postPriorOverlap(tst, dbeta, shape1=0.2, shape2=0.2)

# check overlap with a Uniform(0, 1) prior:
postPriorOverlap(tst, runif(1e6))

}
