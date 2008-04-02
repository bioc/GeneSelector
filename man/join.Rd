\name{join}
\alias{join}
\title{Combine two objects of class RepeatRankings}
\description{
  Convenience method to combine several objects of class
  \link{RepeatRanking}, typically with different Resampling
  types, but based on the same data. The output is again
  an object of class \code{RepeatRanking} (repeated application
  is therefore possible). Useful if one is not exclusively
  in one resampling scheme.
}
\usage{
join(RR1, RR2)
}

\arguments{
  \item{RR1,RR2}{Objects of class \code{RepeatRanking}.}
}

\value{An object of class \code{RepeatRanking}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{RepeatRanking}}
\keyword{univar}
