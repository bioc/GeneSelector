\name{GetStabilityGLM}
\alias{GetStabilityGLM}
\title{Stability measures for significance findings}
\description{
  Assesses the stability of the set of genes declared statistically
  significant for differential expression. To this end, p-values
  or adjusted p-values are used to generate binary response variables
  for a logistic regression model. As single covariate, the ranks
  obtained from the original dataset are used.
  Analogously to the linear model approach, weights are incorporated
  to attribute more importance to higher ranked genes.
  The deviance(s) resulting from these models are used as stability measure.
}
\usage{
GetStabilityGLM(RR, decay = c("linear", "quadratic", "exponential"),
                scheme = c("rank", "pval"), alpha = 1,
                maxpval = 0.05, method=c("raw", "BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{An object of class \link{RepeatRanking}.}
  \item{scheme}{Whether ranks (\code{scheme="rank"}) or p-values (\code{scheme="pval"})
                should be used as basis of weighting.}
  \item{decay}{argument controlling the weight decay for the weights 
               used in the linear regression model. If \code{decay=linear},
               then the weight of the s-th rank/p-value is \code{1/s},
               if \code{decay=quadratic}, then the weight is \code{1/s^2}
               and if \code{decay=exponential}, then the weight is \code{exp(-s*alpha)},
               where \code{alpha} is a tuning parameter specified via
               the argument \code{alpha}.}
  \item{alpha}{To be specified only if \code{decay="exponential"}, s. also
               \link{GetAlpha}.}
  \item{maxpval}{The maximum p-value that is still considered significant
                 (type I error). Default is 0.05.}
   \item{method}{The method used for p-value adjustment, s. \link{AdjustPvalues}.
                If \code{method = "raw"}, then the raw p-values will be used.}
}

\value{An object of class \code{GetStabilityGLM}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, \link{GetStabilityOverlap},
          \link{GetStabilityPCA}, \link{RecoveryScore}, \link{GetAlpha}}
\keyword{univar}

