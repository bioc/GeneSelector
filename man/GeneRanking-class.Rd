\name{GeneRanking-class}
\docType{class}
\alias{GeneRanking-class}
\alias{GeneRanking}
\alias{show,GeneRanking-method}
\alias{summary,GeneRanking-method}
\alias{toplist,GeneRanking-method}

\title{"GeneRanking"}
\description{Object returned by all implemented Ranking methods (\link{RankingTstat},\link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}, \link{RankingGap}) 
             }

\section{Slots}{
	 \describe{
    \item{\code{x}:}{A \code{matrix} storing gene expression,
                     rows correspond to genes, columns to samples (arrays).}
    \item{\code{y}:}{A \code{factor} with two levels of class labels.}
    \item{\code{statistic}:}{A \code{numeric} vector storing the statistics,
                             ordered according to the \code{ranking}.}
    \item{\code{ranking}:}{A vector of indices that represents the
                           ranking of the genes. The first entry corresponds
                           to the best one. The ranking is determined
                           via the statistic (size of absolute value), so
                           there is no distinction of over- and underexpression.
                           For example, if there are five genes with indices
                           \code{1,2,3,4,5}, \code{ranking=3,4,2,1,5} means
                           that the gene with index 3 is ranked highest, having
                           the largest statisitic (in absolze value), gene 4
                           has rank 2, and so on.}
    \item{\code{pval}:}{The ordered vector of p-values.\cr
                        \code{NA} if pvalues have not been computed.}
    \item{\code{type}:}{Type of the test (one of "\code{unpaired}", "\code{paired}",
                         "\code{onesample}").}
    \item{\code{method}:}{Short name of the ranking Method.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(GeneRanking-object)} for brief information.}
    \item{\link[summary]{summary,GeneRanking}}{Use \code{summary(GeneRanking-object)}
                                            For a five-point-summary of statistics and p-values.}
          
    \item{toplist}{Use \code{toplist(object, k=10)} 
                    to get information about the top \code{k=10} genes.}
	 }
}

\keyword{univar}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{RepeatRanking}, \link{GetRepeatRanking}}
