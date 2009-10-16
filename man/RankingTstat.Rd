\name{RankingTstat}
\alias{RankingTstat}
\alias{RankingTstat-methods}
\alias{RankingTstat,matrix,numeric-method}
\alias{RankingTstat,matrix,factor-method}
\alias{RankingTstat,ExpressionSet,character-method}
\title{Ranking based on the 'ordinary' t statistic.}
\description{
  Performs univariate (rowwise) t tests on a gene expression matrix.
}
\usage{
RankingTstat(x, y, type = c("unpaired", "paired", "onesample"), pvalues = TRUE, gene.names = NULL, ...)
}
\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.\cr
           If \code{type = "paired"}, the first half of the columns corresponds to
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
           If \code{type = "paired"}, take care that the coding is
           analogously to the requirement concerning \code{x}.}.
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test, equal variances assumed. \cr
                                   For unequal variances, use \link{RankingWelchT}.}
                \item{"paired":}{paired test. Take care that the coding of \code{y}
                                 is correct (s. above).}
                \item{"onesample":}{\code{y} has only one level. 
                                    Test whether the true mean is different
                                    from zero.}
                   
                                   
                }
                }
  \item{pvalues}{Should p-values be computed ? Default is \code{TRUE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}

\value{An object of class \link{GeneRanking}.}
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\seealso{
         \link{RepeatRanking}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingTstat
ordT <- RankingTstat(xx, yy, type="unpaired")
}
