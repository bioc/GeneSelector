\name{RecoveryScore}
\alias{RecoveryScore}
\title{Stability measures for gene rankings}
\description{
 Computes the Recovery Score of Pavlidis et al. (s. reference)
 below. The stability meausre is the proportion of genes that are 
 declared significant (using usually multiple testing procedures) 
 in both the original and the perturbed dataset.
}
\usage{
RecoveryScore(RR, method = c("raw", "BH", "qvalue", "Bonferroni", "Holm", 
              "Hochberg", "SidakSS", "SidakSD", "BY"), maxpval = 0.05)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}.}
  \item{method}{The p-value adjustment method, s. \link{AdjustPvalues}.
                Can also be \code{"raw"}(default), then no adjustment
                will be done.}
  \item{maxpval}{The maximum p-value at which a gene is still considered
                 significantly differentially expressed (after adjustment).}
}

\value{A numeric vector of recovery scores for each perturbed dataset.}
\references{Pavlidis, P., Li, Q., Noble, W.S. (2003). \cr
            The effect of replication on gene expression microarray experiments.
            \emph{Bioinformatics, 19, 1620-1627}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GetStabilityLm}, \link{GetStabilityOverlap}}
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
### Compute Recovery Score
rs_ordT <- RecoveryScore(loor_ordT, method="BH")
}
