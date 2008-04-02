\name{RankingPermutation}
\alias{RankingPermutation}
\title{Ranking based on permutation tests.}
\description{
 The function is a wrapper for \code{mt.sample.teststat}
 from the package \code{multtest} (Dudoit et al, 2003). The ranking is here
 based on permutation p-values first, followed by the
 absolute value of the statistic.
 For \code{S4} method information, see \link{RankingPermutation-methods}.   
}
\usage{
RankingPermutation(x, y, type = "unpaired", B = 100, gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.}
  \item{type}{Only the two sample case, \code{type="unpaired"} is possible.}
  \item{B}{The number of permutations to generate. Defaults to 100,
           but should be increased if computing power admits. Taking
           \code{B} too high, however, can lead to long computation
           time, especially if called from \link{GetRepeatRanking}}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Further arguments passed to \code{mt.sample.teststat}
               from the package \code{multtest}. Can be used, for example,
               to select the statistic to be computed. By default
               this is \code{"t.equalvar"} (t-test with equal variances assumed).}
}

\note{The p-values, on which the ranking is primarily based, suffer from
      the discreteness of the procedure. They follow a step function
      with jump heights \code{1/B}.}

\value{An object of class \code{GeneRanking}}
\references{Dudoit, S., Shaffer, J.P., Boldrick, J.C. (2003). \cr 
            Multiple Hypothesis Testing in Microarray Experiments
            \emph{Statistical Science, 18, 71-103}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{
         \link{GetRepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
        \link{RankingGap}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingPermutation (100 permutations)
perm <- RankingPermutation(xx, yy, B=100, type="unpaired")
}

