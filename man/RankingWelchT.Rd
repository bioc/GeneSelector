\name{RankingWelchT}
\alias{RankingWelchT}
\alias{RankingWelchT-methods}
\alias{RankingWelchT,matrix,numeric-method}
\alias{RankingWelchT,matrix,factor-method}
\alias{RankingWelchT,ExpressionSet,character-method}
\title{Ranking based on the Welch t statistic.}
\description{
Performs univariate (rowwise) Welch tests on a gene expression matrix.
  The Welch t statistic is a better alternative to the 
  'ordinary' t statistic in the two sample, unequal variances
  setting.
}
\usage{
RankingWelchT(x, y, type = "unpaired", pvalues = TRUE, gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.}
  \item{type}{Only the two sample case, \code{type="unpaired"} is possible.
              Otherwise, use \link{RankingTstat}. Variances are assumed
              to be unequal.}
  \item{pvalues}{Should p-values be computed ? Default is \code{TRUE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currenly unused argument.}
}

\value{An object of class \link{GeneRanking}.}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\seealso{
         \link{RepeatRanking}, \link{RankingTstat}, \link{RankingFC},  \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingWelch
welchT <- RankingWelchT(xx, yy, type="unpaired")
}
