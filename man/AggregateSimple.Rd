\name{AggregateSimple}
\alias{AggregateSimple}
\title{Simple aggregation of repeated rankings}
\description{
  All rankings obtained from perturbed datasets plus the ranking
  from the original dataset are aggregated via \code{aggregatefun}.  
}
\usage{
AggregateSimple(RR, S, aggregatefun = c("mode", "mean", 
                 "median", "quantile"), q = NULL)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}.}
  \item{S}{Either an object of class \code{StabilityLm} or \code{StabilityOverlap}.}
  \item{aggregatefun}{The statistic to return as aggregation.
                      \describe{
                      \item{mode}{The rank occuring most frequently.
                                  If two or more ranks occur equally
                                  often, then weights are used (s.details)}
                      \item{mean}{A weighted mean is used. For information
                                  on weights, s.details.} 
                      \item{median}{The median of all observed ranks is used.}
                      \item{quantile}{The \code{q}-quantile of all observed
                                      ranks is used.}
                      }}
  \item{q}{Only specified if \code{aggregatefun=quantile}}
}
\details{
  The weights used if \code{aggregatefun=mode} or 
  \code{aggregatefun=mean} decrease linear with
  decreasing similarity of perturbed dataset and original dataset
  as measured by Stability Measures (object \code{S}).
}

\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, \link{GetStabilityOverlap},
         \link{AggregateBayes}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingTstat
ordT <- RankingTstat(xx, yy, type="unpaired")
### Generate Leave-one-out Foldmatrix
loo <- GenerateFoldMatrix(xx, yy, k=1)
### Get all rankings
loor_ordT <- GetRepeatRanking(ordT, loo)
### compute stability measure
stab_overlap <- GetStabilityOverlap(loor_ordT, decay="linear")
### aggregate rankings
agg_simple_ordT <- AggregateSimple(loor_ordT, stab_overlap, aggregatefun="mean")
}
