\name{RankingGap}
\alias{RankingGap}
\title{Ranking based on 'gaps'}
\description{
  The ranking is based on the \emph{gap} between two classes
  (for \code{type="unpaired"}) where the gap is defined as
  \code{gap(1,2) = max(min2 - max1, min1-max2, 0)}, where
  min1, max1 are the minimum/maximum observed values from class 1
  (analogously for min2, max2). It is only greater than zero
  if classes do not overlap. For \code{type="paired", "onesample"}
  the gap, i.e. absolute distance from the origin is computed.\cr
  For \code{S4} method information, see \link{RankingFC-methods}.  
}
\usage{
RankingGap(x, y, type = c("unpaired", "paired", "onesample"), gene.names = NULL, ...)
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
                \item{"unpaired":}{two-sample test.}
                \item{"paired":}{paired test. Take care that the coding of \code{y}
                                 is correct (s. above)}
                \item{"onesample":}{\code{y} has only one level. 
                                    Test whether the true mean is different
                                    from zero.}
                }
                }
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}

\note{In most cases, classes will not be separated by only one gene.
      Consequently, the great majority of statistics will be zero.\cr
      p-values are \emph{not} available.}

\value{An object of class \link{GeneRanking}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}


\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingGap
gapstat <- RankingGap(xx, yy, type="unpaired")
}
\keyword{univar}
