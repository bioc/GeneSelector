\name{CombinedRanking-class}
\docType{class}
\alias{CombinedRanking-class}
\alias{CombinedRanking}
\alias{show,CombinedRanking-method}
\alias{toplist,CombinedRanking-method}
\alias{SelectedGenes-methods}
\alias{SelectedGenes}
\alias{SelectedGenes,CombinedRanking-method}
\alias{GeneInfoScreen,CombinedRanking-method}
\alias{GeneInfoScreen-methods}
\alias{GeneInfoScreen}
\title{"CombinedRanking"}

\description{Object returned from a call to \code{GeneSelector}.}

\section{Slots}{
	 \describe{
    \item{\code{ranking}:}{'Final' ranking, usually obtained
                           by combinining several statistics.
                           The generation of this ranking
                           is described in \link{GeneSelector}.\cr
                           The first entry contains the index of
                           the gene ranked highest, as in \link{GeneRanking}.}
    \item{\code{rankmatrix}:}{Matrix of rankings, arranged in a way
                              that the rankings from the most important
                              statistic is in column 1. In contrast to
                              \code{ranking}, the first row contains
                              the rankings of the gene with the first \emph{index}.
                              This slot is used rather for internal reasons.}
    \item{\code{inout}:}{Matrix arranged in the same way as \code{rankmatrix},
                         but information is now binary: If the specified threshold,
                         then there is a \code{"+"} symbolizing selection,
                         whereas a \code{"-"} symbolizes removal.}
    \item{\code{selected}:}{The indices of those genes that fall below
                             the specified threshold. Can be accessed
                             more conveniently using \code{SelectedGenes}}
    \item{\code{adjpval}:}{Numeric vector of adjusted p-values, ordered
                           according to \code{ranking}. \code{NA} if
                           no adjustment has been taken place (in the case
                           that the threshold was fixed by the user).}
    \item{\code{maxrank}:}{Threshold rank, either defined by the user or
                           obtained via p-value adjustment.}
    \item{\code{statistics}:}{The names of the statistics used, ordered
                              according to their importance 
                              (as defined by the user).}
    \item{\code{absdist}:}{Absolute (\code{L1}) distance from (theoretically)
                            best possible result (rank 1 in all rankings),
                            ordered according to \code{ranking}. Note
                            that minimum \code{absdist} does not imply
                            best rank and vice versa, because the computed
                            distance does not weigh different statistics
                            differently.}
    \item{\code{reldist}:}{A 'normalized' version of \code{absdist}
                           (lies in [0;1]).}
  }
}
\section{Methods}{
  \describe{    
    \item{show}{Use \code{show(object)} for brief information.}
    \item{toplist}{Use \code{toplist(object, k=10)}  
          to get information about the top \code{k=10} genes.\cr
          The ranking used is a synthesis from several statistics.}     
    \item{SelectedGenes}{Use \code{SelectedGenes(object)} to show all genes that
          have been selected by the \link{GeneSelector}.}
    \item{GeneInfoScreen}{Use \code{GeneInfoScreen, which=1} to get detailed
          information about the gene with index 1, arranged in
          a pretty plot.}
    \item{plot}{Use \code{plot(object)} to visualize relative distances,
                s. \link{plot,CombinedRanking}}
   }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}

