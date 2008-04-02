\name{AggregatePCA}
\alias{AggregatePCA}
\title{aggregation of repeated rankings by principal components}
\description{
 A principal components analysis is applied to the matrix storing
 the different rankings for each gene. The first principal component
 is then used for aggregation.
}
\usage{
AggregatePCA(RR)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}.}
}


\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\link{AggregateSimple}, \link{AggregateBayes}, \link{AggregatePenalty},
         \link{AggregatePCA}}
\keyword{univar}
