\name{plot,RepeatRanking}
\alias{plot,RepeatRanking}
\alias{plot,RepeatRanking,missing-method}
\title{Visualize results from GetRepeatRanking}
\description{
 A scatterplot of rankings in perturbed datasets (y-axis)
 vs. the ranking of the original dataset (x-axis).
}

%\usage{plot(RepeatRanking-object)}
\arguments{
   \item{x}{An object of class \code{RepeatRanking}.}
   \item{frac}{The fraction of top genes for which the plot is done.
               Default is 1/100.}
   \item{...}{Additional arguments concerning graphical options.}
 }


\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, \link{GetStabilityOverlap}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}