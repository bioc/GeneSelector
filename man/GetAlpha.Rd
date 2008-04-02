\name{GetAlpha}
\alias{GetAlpha}
\title{Helper function for stability assessement.}
\description{
  Both \link{GetStabilityLm} and \link{GetStabilityOverlap}
  depend on a parameter \code{alpha} if \code{decay=exponential}.
  If the weights are based on ranks, then a nonlinear regression
  of the form \cr
  \code{pval = 1- exp(-alpha*rank)} \cr
  can be used to find an appropriate value for \code{alpha}
  via nonlinear least squares. In order to adjust for too
  'optimistic' p-values, multiple testing adjustments should
  be used, s. \link{AdjustPvalues}. 
}
\usage{
GetAlpha(ranking, pval, alpha0 = 0.01)
}
\arguments{
  \item{ranking}{A numeric vector of ranks, regarded as regressor.}
  \item{pval}{A numeric vector of p-values corresponding to the vector
              \code{ranking}.} 
  \item{alpha0}{A starting value for the nonlinear least squares
               estimation procedure passed to \link[=nls]{nls}}
}

\value{The nonlinear least squares estimator for \code{alpha}, s.
       description.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\note{It is more or less  equivalent to use a p-value based ranking
      instead of ranks combined with this procedure.}
\seealso{\link{GetStabilityLm}, \link{GetStabilityOverlap}, \link{AdjustPvalues},
         \link[=nls]{nls}}
\examples{
### rankings
ranks <- 1:100
### corresponding p-values
pvals <- 1-exp(-0.01*ranks) + rnorm(100, sd=0.001)
### determine alpha
alphaopt <- GetAlpha(ranks, pvals, alpha0 = 0.01)
}
\keyword{univar}
