\name{RepeatRanking-class}
\docType{class}
\alias{RepeatRanking-class}
\alias{RepeatRanking}
\alias{show,RepeatRanking-method}
\alias{toplist,RepeatRanking-method}

\title{"RepeatRanking" }
\description{Object returned by a call to \link{GetRepeatRanking}}
\section{Slots}{
	 \describe{
    \item{\code{original}:}{The ranking based on the original data set,
                            represented by an object of class \code{"GeneRanking"}}
    \item{\code{rankings}:}{The rankings obtained from altered datasets,
                            stored as a matrix.
                            One column represents one replication. Each
                            column is arranged in the same manner as the
                            slot \code{ranking} of the class \code{GeneRanking},
                            i.e. the first entry of a column contains the
                            index of the gene ranked highest.}
    \item{\code{pvals}:}{The p-values obtained from altered data sets,
                         stored analogoulsy to \code{rankings}. If p-values
                         have not been computed, this is a matrix of \code{NA}s.}
    \item{\code{statistics}:}{The statistics obtained from altered data sets,
                              stored analogoulsy to \code{rankings}}
    \item{\code{scheme}:}{A character for the resampling scheme, can be one
                          of \code{"Subsampling"}, \code{"Labelexchange"},
                          \code{"Bootstrap"}, \code{"Jittering"} (if noise
                          has been added) or \code{"combined"} if several
                          resampling schemes for the same dataset and ranking method
                          have been combined via the \code{join}-method, s. below.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{use \code{show(RepeatRanking-Object)} for brief information.}
    \item{toplist}{Use \code{toplist(RepeatRanking-Object, k=10)}  
                   to get information about the top \code{k=10} genes
                   for each replication (=perturbed dataset) and one
                   overall table showing frequencies of gene indices
                   for each of the ranks \code{1,...k}. 
                   Additionally, only the overall table can be shown with
                   all other output suppressed using 
                   \code{toplist(RepeatRanking-Object, show=FALSE)}}
    \item{variance}{Genewise variance estimation, s. \link{variance,RepeatRanking-method}} 
    \item{join}{use \code{join(RepeatRanking-Object1, RepeatRanking-Object2)}
                to combine results from different resampling schemes. The
                results is again an object of class \code{RepeatRanking}
                where the slot \code{scheme} is \code{"combined"} and
                all matrices have been concatenated columnwise.}
    \item{plot}{use \code{plot(RepeatRanking-Object)} for a scatterplot
                of original rankings and rankings of the perturbed datasets,
                s. \link{plot,RepeatRanking}}
   
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GeneRanking}, \link{GetRepeatRanking}}
\keyword{univar}