\name{RankingWilcEbam}
\alias{RankingWilcEbam}
\title{Ranking based on the empirical bayes approach of Efron}
\description{
  The function is a wrapper for the function \code{wilc.ebam}
  from the package \code{siggenes} that implements an
  empirical bayes mixture model approach in combination
  with the Wilcoxon statistic.\cr
   For \code{S4} method information, see \link{RankingWilcEbam-methods}.
}
\usage{
RankingWilcEbam(x, y, type = c("unpaired", "paired", "onesample"), gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet} alternatively an object of class \code{ExpressionSet}.\cr
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
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Further arguments to be passed to \code{wilc.ebam},
         s. package \code{siggenes}.}
}

\note{p-values are \emph{not} computed - the statistic is
      a posterior probabiliy.}

\value{An object of class \code{GeneRanking}.}
\references{Efron, B., Tibshirani, R. (2002). \cr
            Empirical Bayes Methods and False Discovery Rates for Microarrays
            \emph{Genetic Epidemiology, 23, 70-86}
            
            Schwender, H., Krause, A. and Ickstadt, K. (2003).\cr
             Comparison of the Empirical Bayes and the Significance 
             Analysis of Microarrays. 
             \emph{Techical Report, University of Dortmund.}
             }

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam},  \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}, \link{RankingGap}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingWilcEbam
WilcEbam <- RankingWilcEbam(xx, yy, type="unpaired")
}
