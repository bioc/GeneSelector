\name{plot,AggregatedRanking}
\alias{plot,AggregatedRanking}
\alias{plot,AggregatedRanking,missing-method}
\title{Visualize results from AggregateBayes}
\description{Display of the (discrete) posterior distribution
             of the rank (with respect to differential expression)
             of a certain gene}

%\usage{plot(x,...)} 
\note{Only works if the aggregation has been done with \code{AggregateBayes}.
For \code{AggregateSimple}, there is no plot method.}
\arguments{\item{x}{An object of class \code{AggregateBayes}}
           \item{index}{the gene index (row of expression matrix)
                        for which is the plot is performed}
           \item{...}{Additional arguments concerning graphical options.}}


\seealso{\link{GetRepeatRanking}, \link{AggregateBayes}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}