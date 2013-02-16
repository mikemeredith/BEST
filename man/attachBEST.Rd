\name{attachBEST}
\alias{attachBEST}
\alias{detachBEST}
\title{
Attach / detach elements of BEST objects to search path
}
\description{
The BEST object is attached/detached to the search path, so that the samples for each parameter can be accessed by name: mu1, mu2, nu, sigma1, sigma2 for the 2-group case, or mu, nu, sigma for a 1-group analysis.
}
\usage{
attachBEST(BESTobj, overwrite = NA)
detachBEST()
}
\arguments{
  \item{BESTobj}{
an object of class 'BEST', as produced by the function \code{\link{BESTmcmc}}.
}
  \item{overwrite}{
logical: should objects in the Global environment with the same names (eg. mu1, nu or sigma) be removed? If NA, the default, a dialogue box will appear asking if objects should be removed.
}
}
\details{
\code{attachBEST} first creates a data frame with one column for each parameter and gives the columns appropriate names. The resulting data frame is then attached to the search path. See \code{\link{attach}} for details.

Only one BEST object can be attached at one time; if \code{attachBEST} is used again, the previous object is detached before the new one is attached.

\code{detachBEST} removes the BEST object (if any) from the search path.
}
\value{
attachBEST returns invisibly the environment with a "name" attribute.

detachBEST returns NULL.
}
\author{
Mike Meredith
}
\seealso{
\code{\link{attach}}. For examples, see \code{\link{BEST-package}}.
}
\examples{
# See the examples in BEST-package.
}
