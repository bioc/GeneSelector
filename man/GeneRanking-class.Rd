\name{GeneRanking-class}
\docType{class}
\alias{GeneRanking-class}
\alias{GeneRanking}
\alias{show,GeneRanking-method}
\alias{summary,GeneRanking-method}
%\alias{summary,GeneRanking}
\alias{toplist,GeneRanking-method}

\title{"GeneRanking"}
\description{Object returned by all implemented ranking methods (\link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}).}

\section{Slots}{
  \describe{
    \item{\code{x}:}{Gene expression\code{matrix} rows correspond to genes, columns to samples (arrays).}
    \item{\code{y}:}{A two-level \code{factor} of class labels.}
    \item{\code{statistic}:}{A \code{numeric} vector storing the test statistics.}
    \item{\code{ranking}:}{The ranking is determined via the statistic (normally via the size of absolute value).
                           The lower the rank, the higher the magnitude of differential expression.
                           There is no distinction of over- and underexpression.}
    \item{\code{pval}:}{The vector of p-values computed from \code{statistic}.\cr
                        \code{NA} if p-values have not been computed.}
    \item{\code{type}:}{Type of the test (one of "\code{unpaired}", "\code{paired}",
                         "\code{onesample}").}
    \item{\code{method}:}{Short name of the ranking method.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(GeneRanking-object)} for brief information.}
    \item{summary}{Use \code{summary(GeneRanking-object)} for a five-point-summary of statistics and p-values arranged as two-column table. The second column (p-values) are \code{NA}
  in the case that p-values have not been computed.}
    \item{toplist}{Use \code{toplist(object, k=10)} 
                    to get information about the top \code{k=10} genes.}
	 }
}

\keyword{univar}
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}

\seealso{\link{RepeatRanking}}
