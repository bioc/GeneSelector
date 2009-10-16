\name{RankingEbam}
\alias{RankingEbam}
\alias{RankingEbam-methods}
\alias{RankingEbam,matrix,numeric-method}
\alias{RankingEbam,matrix,factor-method}
\alias{RankingEbam,ExpressionSet,character-method}
\title{Ranking based on the empirical Bayes approach of Efron et al. (2001)}
\description{
  The approach of Efron et al. (2001) is based on a mixture model
  for two subpopulations: genes that are differentially expressed and
  those not. The posterior probability for differential
  expression is used to obtain a ranking.
  The function described below is merely a wrapper for the
  function \code{z.ebam} from the package \code{siggenes}.\cr
  For \code{S4} method information, see \link{RankingEbam-methods}.  
}
\usage{
RankingEbam(x, y, type = c("unpaired", "paired", "onesample"), gene.names = NULL, ...)
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
           analogously to the requirement concerning \code{x}.}
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test.}
                \item{"paired":}{paired test. Take care that the coding of \code{y}
                                 is correct (s. above)}
                \item{"onesample":}{\code{y} has only one level. 
                                    Test whether the true mean is different
                                    from zero.}
                   
                                   
                }}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Further arguments passed to the function \code{z.ebam}.
               Can be used to influence the \emph{fudge factor} to
               the stabilize the variance. Currently, the 90 percent
               quantile is used.}
}
\details{
  To find a better value for the fudge factor, the function 
  \code{find.a0} (package \code{siggenes}) can be used.
}
\note{P-values are \emph{not} computed - the statistic is
      a posterior probabiliy.}
\value{An object of class \link{GeneRanking}.}
\references{Efron, B., Tibshirani, R., Storey, J.D., Tusher, V. (2001).\cr 
             Empirical Bayes Analysis of a Microarray Experiment. 
             \emph{Journal of the American Statistical  Association, 96, 1151-1160}.
              
             Schwender, H., Krause, A. and Ickstadt, K. (2003).\cr
             Comparison of the Empirical Bayes and the Significance 
             Analysis of Microarrays. 
             \emph{Techical Report, University of Dortmund.}}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\seealso{
         \link{RepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
      \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingEbam
Ebam <- RankingEbam(xx, yy, type="unpaired")
}
