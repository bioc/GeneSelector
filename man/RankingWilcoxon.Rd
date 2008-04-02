\name{RankingWilcoxon}
\alias{RankingWilcoxon}
\title{Ranking based on the Wilcoxon statistic}
\description{
  The Wilcoxon statistic is rank-based and 'distribution free'. 
It is equivalent to the Mann-Whitney statistic and also related to the 'Area under
  the curve' (AUC) in the two sample case. The implementation
  is efficient, but still far slower than that of the t-statistic.\cr
  For \code{S4} method information, see \link{RankingWilcoxon-methods}.  
}
\usage{
RankingWilcoxon(x, y, type = c("unpaired", "paired", "onesample"), pvalues = FALSE, gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.\cr
           If \code{type = paired}, the first half of the columns corresponds to 
           the first measurements and the second half to the second ones. 
           For instance, if there are 10 observations, each measured twice,
           stored in an expression matrix \code{expr}, 
           then \code{expr[,1]} is paired with \code{expr[,11]}, \code{expr[,2]}
           with \code{expr[,12]}, and so on.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.\cr
           If \code{type = paired}, take care that the coding is
           analogously to the requirement concerning \code{x}}
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test, Wilcoxon Rank Sum test
                                   is performed.}
                \item{"paired":}{Wilcoxon sign rank test is performed
                                 on the differences.}
                \item{"onesample":}{\code{y} has only one level. 
                                    The Wilcxon sign rank test for difference from zero is
                                    performed.}
                
                 }                  
                }
  \item{pvalues}{Should p-values be computed ? Default is \code{FALSE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}

\note{Note that although the Wilcoxon Rank Sum test is distribution-free,
      it is not without assumptions.}

\value{An object of class \link{GeneRanking}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, 
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}, \link{RankingGap}, \link[=wilcox.test]{wilcox.test}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingWilcoxon
wilcox <- RankingWilcoxon(xx, yy, type="unpaired")
}

