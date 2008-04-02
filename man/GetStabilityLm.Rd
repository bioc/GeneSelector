\name{GetStabilityLm}
\alias{GetStabilityLm}
\title{Stability measures for gene rankings}
\description{
  Assesses stability of gene rankings by regressing
  the rankings of perturbed datasets on the ranking
  of the original datasets in a weighted manner.
  The idea is that if stability is high, the resulting
  regression models fit well.
}
\usage{
GetStabilityLm(RR, decay = c("linear", "quadratic", "exponential"), 
                measure = c("wilks", "direct"), 
                scheme = c("rank", "pval"), alpha = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{An object of class \link{RepeatRanking}.}
  \item{scheme}{Whether ranks (\code{scheme="rank"})
                or p-values (\code{scheme="pval"})
                should be used as basis of weighting.}
  \item{decay}{argument controlling the weight decay for the weights 
               used in the linear regression model. If \code{decay=linear},
               then the weight of the s-th rank/p-value is \code{1/s},
               if \code{decay=quadratic}, then the weight is \code{1/s^2}
               and if \code{decay=exponential}, then the weight is \code{exp(-s*alpha)},
               where \code{alpha} is a tuning parameter specified via
               the argument \code{alpha}.}
  \item{measure}{The stability measure to be computed. If \code{measure="wilks"},
                 then a stability measure based on the Wilk's Lambda Test
                 for multivariate linear regression models is used. If 
		\code{measure="direct"}, then the direct generalization
                 of the univariate coefficient of determination to the multivariate
                 case is used. The second approach can fail if there exists rankings of the
                 perturbed datasets that are exactly equal.}
  
  \item{alpha}{To be specified only if \code{decay="exponential"}, s. also
               \link{GetAlpha}.}
  \item{\dots}{Further arguments passed to \code{lm}.}
}

\value{An object of class \code{GetStabilityLm}}
\references{Mardia, K.V., Kent, J.T., Bibby, J.M. (1979). \cr
            Multivariate Analysis \emph{Academic Press.}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GetRepeatRanking}, \link{GetStabilityOverlap}, 
          \link{RecoveryScore}, \link{GetAlpha}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### get ranking 
ordT <- RankingTstat(xx, yy, type="unpaired")
### Generate Leave-One-Out
loo <- GenerateFoldMatrix(xx, yy, k=1)
### Repeat Ranking with t-statistic
loor_ordT <- GetRepeatRanking(ordT, loo)
### assess stability
stab_lm_ordT <- GetStabilityLm(loor_ordT, decay="linear")
### plot
plot(stab_lm_ordT )
}
