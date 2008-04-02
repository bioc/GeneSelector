\name{plot,CombinedRanking}
\alias{plot,CombinedRanking}
\alias{plot,CombinedRanking,missing-method}
\title{Visualize results from GeneSelector}
\description{The bars in this barplot symbolize the L1 (absolute)
             distance from the best possible results
             (rank 1 for all statistics). The ordering on
             the axis can disagree with the heights of the bars
             due to the fact that all statistics are equally
             weighted independent of the order of different
             statistics defined for the call to \code{GeneSelector}.}

%\usage{plot(CombinedRanking-object)}
\arguments{\item{x}{An object of class \code{CombinedRanking}.}
           \item{top}{the top number of genes for which the plot is done.}
           \item{...}{Additional arguments concerning graphical options.}}


\seealso{\link{GeneSelector}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}

