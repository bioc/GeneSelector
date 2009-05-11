\name{RepeatedRanking-class}
\docType{class}
\alias{RepeatedRanking-class}
\alias{RepeatedRanking}
\alias{show,RepeatedRanking-method}
\alias{toplist,RepeatedRanking-method}
\alias{plot,RepeatedRanking}
\alias{plot,RepeatedRanking,missing-method}
\alias{MergeRankings}
\alias{MergeRankings-methods}
\alias{MergeRankings,RepeatedRanking,RepeatedRanking-method}

\title{"RepeatedRanking"}
\description{Object returned by a call to \link{RepeatRanking} or \link{MergeMethods}}
\section{Slots}{
	 \describe{
    \item{\code{original}:}{The ranking based on the original data set (output from \link{RepeatRanking}) or on a reference
                            method (output from \link{MergeMethods}), represented by an object of class \code{GeneRanking}}
    \item{\code{rankings}:}{The rankings obtained from altered datasets (output from \link{RepeatRanking} or from different methods (output from \link{MergeMethods}),
                            stored as a matrix.
                            One column represents one replication (output from \link{RepeatRanking} or one method (output from \link{MergeMethods}). Each
                            column is of the same structure as the
                            slot \code{ranking} of the class \code{GeneRanking}.}
    \item{\code{pvals}:}{The matrix of p-values stored analogously to \code{rankings}. If p-values
                         have not been computed, this is a matrix of \code{NA}s.}
    \item{\code{statistics}:}{The statistics obtained from altered data sets,
                              stored analogously to \code{rankings}}
    \item{\code{scheme}:}{A character for the resampling scheme, can be one
                          of \code{"subsampling"}, \code{"labelexchange"},
                          \code{"bootstrap"}, \code{"jittering"} (if noise
                          has been added), \code{"merged (rankings)"} if several
                          resampling schemes for the same dataset and ranking method
                          have been combined via the \code{MergeRankings}-method, or
                          \code{"merged (methods)"} if the rankings of several methods have been
                          combined via the \code{MergeMethods}-method.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(RepeatedRanking-Object)} for brief information.}
    \item{toplist}{Use \code{toplist(RepeatedRanking-Object, k=10)}
                   to get information about the top \code{k=10} genes
                   for each ranking and one
                   overview table showing frequencies of gene indices
                   for each of the ranks \code{1,...k}. 
                   Additionally, only the overview table can be shown with
                   all other output suppressed using 
                   \code{toplist(RepeatedRanking-Object, show=FALSE)}.}
    \item{dispersion}{Genewise variance estimation, s. \link{dispersion,RepeatedRanking-method}}
    \item{MergeRankings}{Use \code{MergeRankings(RepeatedRanking-Object1, RepeatedRanking-Object2)}
                to combine results from different resampling schemes. The
                results is again an object of class \code{RepeatedRanking}
                where the slot \code{scheme} is \code{"merged (rankings)"} and
                all matrices have been concatenated columnwise.}
    \item{plot}{Use \code{plot(RepeatedRanking-Object)} for a scatterplot
                of the reference ranking (slot \code{original}) vs. alternative rankings (slot \code{rankings}).}
    \item{HeatmapRankings}{s. \link{HeatmapRankings}}
   
	 }
}

\author{Martin Slawski \cr
        Anne-Laure Boulesteix}

\seealso{\link{GeneRanking}, \link{RepeatRanking}, \link{MergeMethods}, \link{dispersion}, \link{HeatmapRankings}}
\keyword{univar}