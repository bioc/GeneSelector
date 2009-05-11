\name{RankingShrinkageT}
\alias{RankingShrinkageT}
\alias{RankingShrinkageT-methods}
\alias{RankingShrinkageT,matrix,numeric-method}
\alias{RankingShrinkageT,matrix,factor-method}
\alias{RankingShrinkageT,ExpressionSet,character-method}
\title{Ranking based on the 'shrinkage t' statistic}
\description{
  The shrinkage t statistic stabilizes the estimated variances
  appearing in the denominator of the statistic via a 
  James-Stein-Shrinkage approach (Opgen-Rhein and Strimmer, 2007). In this implementation,
  the shrinkage target is the median of the variances.}
\usage{
RankingShrinkageT(x, y, type = c("unpaired", "paired", "onesample"), gene.names = NULL, ...)
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
           If \code{type = "paired"}, take care that the coding is
           analogously to the requirement concerning \code{x}.}
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test.}
                \item{"paired":}{paired test. Take care that the coding of \code{y}
                                 is correct (s. above).}
                \item{"onesample":}{\code{y} has only one level. 
                                    Test whether the true mean is different
                                    from zero.}
                   
                                   
                }}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}

\value{An object of class \link{GeneRanking}.}
\references{Opgen-Rhein, R., Strimmer, K. (2007). \cr
            Accurate Ranking of Differentially Expressed Genes 
            by a Distribution-Free Shrinkage Approach.
            \emph{Statistical Applications in Genetics and
                  Molecular Biology, Vol. 6, Iss. 1, Art.9}}

\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}

\seealso{
         \link{RepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingShrinkageT
shrinkaget <- RankingShrinkageT(xx, yy, type="unpaired")
}
