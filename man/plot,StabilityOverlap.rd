\name{plot,StabilityOverlap}
\alias{plot,StabilityOverlap}
\alias{plot,StabilityOverlap,missing-method}
\title{Visualize results fromGetStabilityOverlap}
\description{Plots cumulated (top) and averaged overlap (bottom)
             score in dependency of ranks. The bold line
             in the top display depicts the maximum possible
             score. The averaged overlap score (bottom) is at most 1
             and at least 0.}

%\usage{plot(StabilityOverlap-object)}
\arguments{
 \item{x}{An object of class \code{StabilityOverlap}} 
 \item{frac}{The fraction of top genes for which the plot is done.
               Default is 1/50.}
 \item{...}{Additional arguments concerning graphical options.}
 }

\seealso{\link{GetRepeatRanking}, \link{GetStabilityOverlap}}
\references{Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}