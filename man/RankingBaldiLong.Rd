\name{RankingBaldiLong}
\alias{RankingBaldiLong}

\title{Ranking based on the t-statistic of Baldi and Long}
\description{
  Performs bayesian t tests on a gene expression matrix.\cr
  For \code{S4} method information, see \link{RankingBaldiLong-methods}.}
\usage{
RankingBaldiLong(x, y, type = c("unpaired", "paired", "onesample"), 
                  m = 100, conf = NULL, pvalues = TRUE, gene.names = NULL, ...)
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
  \item{m}{Size of the sliding window that is used obtain the background variance
           from pooled similarly expressed genes. s. Details.}
  \item{conf}{The number of 'pseudocounts' giving weight to the prior variance. s. Details.}
  \item{pvalues}{Should p-values be computed ? Default is \code{TRUE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}
\details{
The argument \code{m} determines the width of the window
used to provides an estimate of the average variability of 
gene expression for those genes that show a similar expression level.\cr
The argument \code{conf} is non-negative and indicates the weight give to the Bayesian prior estimate of within-treatment 
variance. Baldi and Long report reasonable performance with this parameter set equal to approximately 3 times the number of 
observations, when the number of experimental observations is small (approximately 4 or less). 
If the number of replicate experimental observations is large then the confidence value can be lowered 
to be equal to the number of observations (or even less).  
}
\value{An object of class \link{GeneRanking}.}
\references{Baldi,P., Long, A.D. (2001). \cr
            A bayesian framework for the analysis of microarray data.
            \emph{Bioinformatics, 17, 509-519}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\note{Results can differ slighlty from the Cyber-T-Software
      of Baldi and Long.}
\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
         \link{RankingFoxDimmic}, \link{RankingLimma}, 
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
### run RankingBaldiLong
BaldiLong <- RankingBaldiLong(xx, yy, type="unpaired")
}
