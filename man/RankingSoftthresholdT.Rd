\name{RankingSoftthresholdT}
\alias{RankingSoftthresholdT}
\alias{RankingSoftthresholdT-methods}
\alias{RankingSoftthresholdT,matrix,numeric-method}
\alias{RankingSoftthresholdT,matrix,factor-method}
\alias{RankingSoftthresholdT,ExpressionSet,character-method}
\title{Ranking via the 'soft-threshold' t-statistic}
\description{
 The 'soft-threshold' statistic (Wu, 2005) is constructed
 using a linear regression model with an \code{L1}
 penalty (also referred to as LASSO penalty). In special
 cases (like here) the LASSO estimator can be 
 calculated analytically and is then known as 'soft threshold
 estimator'.}
\usage{
RankingSoftthresholdT(x, y, type = c("unpaired", "paired", "onesample"),
                      lambda = c("lowess", "cor", "user"), userlambda = NULL, 
                      gene.names = NULL, ...)
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
  \item{lambda}{s. details}
  \item{userlambda}{A user-specified value for \code{lambda}, s. details.}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Currently unused argument.}
}
\details{
  There are currently three ways of specifying the shrinkage intensity
  \code{lambda}. Both \code{"lowess"} and \code{"cor"} are relatively
  slow, especially if rankings are calculated repeatedly (\link{RepeatRanking}).
  Therefore, a 'reasonable' value can be set by the user.
}

\value{An object of class \code{GeneRanking}.}
\references{Wu, B. (2005). Differential gene expression using
            penalized linear regression models: The improved
            SAM statistic. \emph{Bioinformatics, 21, 1565-1571}}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\note{The code is a modified version of an implementation available in the \code{st}
      package of Opgen-Rhein and Strimmer (2007).}
\seealso{
         \link{RepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT},  
          \link{RankingPermutation}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingSoftthresholdT
softt <- RankingSoftthresholdT(xx, yy, type="unpaired")
}



