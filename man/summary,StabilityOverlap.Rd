\name{summary,StabilityOverlap}
\alias{summary,StabilityOverlap}

\title{Summarize Overlap Scores}

\description{Returns a five-point summary of overlap counts
             and scores, where the summary is done with respect
             to the different perturbed datasets.}

%\usage{summary(StabilityOverlap-object)}
\arguments{
 \item{object}{An object of class \code{StabilityOverlap}}
 \item{which}{\describe{
              \item{"scores"}{Summary with respect to  overlap score.}
              \item{"overlap"}{Summary with respect to the number of overlaps.}
              }}
 \item{position}{Up to which position in the list (rank) overlap and score
                 are computed. For instance, if \code{position=100}(default) then
                 the summary is based on overlaps/scores on the lists
                 of the top 100 genes.}
 }


\seealso{\link{StabilityOverlap}}
\references{Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}
