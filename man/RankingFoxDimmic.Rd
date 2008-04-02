\name{RankingFoxDimmic}
\alias{RankingFoxDimmic}
\title{Ranking based on the t-statistic of Fox and Dimmic}
\description{
  Performs a two-sample bayesian t test on a gene expression matrix using
the methodology by Fox and Dimmic (2006).\cr
  For \code{S4} method information, see \link{RankingFoxDimmic-methods}.}
\usage{
RankingFoxDimmic(x, y, type = "unpaired", m = 8, pvalues = TRUE, gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.}
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test, equal variances assumed.}
                }
                \code{"paired"} and \code{"unpaired"} are not possible
                for this kind of test.}
  \item{m}{The number of similarly expressed genes to use for calculating Bayesian variance
           and prior degrees of freedom. The default value suggested
           by Fox and Dimmic is currently 8, s. note.}
  \item{pvalues}{Should p-values be computed ? Default is \code{TRUE}.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}

\value{An object of class \link{GeneRanking}.}
\references{Fox, R.J., Dimmic, M.W. (2006). \cr
            A two sample Bayesian t-test for microarray data.
            \emph{BMC Bioinformatics, 7:126}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\note{Although the test of Fox and Dimmic is very similar to the
     one proposed by Baldi and Long, there are various slight differences,
     in particular with respect to the computation of the bayesian variance.}
\seealso{
        \link{GetRepeatRanking},  \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong},  \link{RankingLimma}, 
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
### run RankingFoxDimmic
FoxDimmic <- RankingFoxDimmic(xx, yy, type="unpaired")
}
