\name{RankingFC}
\alias{RankingFC}
\title{Ranking based on the (log) foldchange}
\description{
  Naive ranking that only considers difference in means
  without taking variances into account.\cr
  For \code{S4} method information, see \link{RankingFC-methods}.  
}
\usage{
RankingFC(x, y, type = c("unpaired", "paired", "onesample"), 
          pvalues = TRUE, gene.names = NULL, LOG = FALSE, ...)
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
                                       
                }}
  \item{pvalues}{Should p-values be computed ? Defaults to \code{TRUE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{LOG}{By default, the data are assumed to be already logarithm-ed. 
             If not, this can be done by setting \code{LOG=TRUE}}
  \item{\dots}{Currently unused argument.}
}

\note{Take care that the \emph{log} foldchange is computed, therefore
logarithmization might be necessary.\cr
The p-values for the difference in means are based on a standard
normal assumption. 
}

\value{An object of class \link{GeneRanking}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat},  \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}, \link{RankingGap}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingFC
FC <- RankingFC(xx, yy, type="unpaired")
}
