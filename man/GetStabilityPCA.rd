\name{GetStabilityPCA}
\alias{GetStabilityPCA}

\title{Stability measures for gene rankings}
\description{
  A principal components analysis is applied to the matrix storing
 the different rankings for each gene. The ratio of the first
 eigenvalue to the sum of all eigenvalues is used as stability measure.
 If stability is high/variability is low, then the first principal
 component will explain a large amount of the overall variation,
 leading to large first eigenvalue.
}

\usage{
GetStabilityPCA(RR)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}}
}
\value{An object of class \link{GetStabilityPCA}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}



\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, \link{GetStabilityOverlap}}
\keyword{univar}