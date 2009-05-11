\name{AggregateSimple}
\alias{AggregateSimple}
\alias{AggregateSimple-methods}
\alias{AggregateSimple,RepeatedRanking-method}
\title{Simple aggregation of repeated rankings}
\description{
  All obtained rankings are aggregated by a genewise summary \code{measure}.
}
\usage{
  AggregateSimple(RR, measure = c("mode", "mean", "trimmed.mean", "median", "quantile"), q=NULL, trim = NULL)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}}
  \item{measure}{The statistic to be used as basis for the aggregated ranking.
                      \describe{
                      \item{mode}{The rank occuring most frequently. If several ranks occur equally often, the lowest one is used.}
                      \item{mean}{The mean of the ranks.}
                      \item{trimmed.mean}{The trimmed mean of the ranks, i.e. the mean resulting when throwing away the \code{trim}*100 percent most extreme observations at both tails.}
                      \item{median}{The median of the ranks.}
                      \item{quantile}{The \code{q}-quantile, \code{0 <= q <= 1}, of the ranks.}}}
  \item{q}{Only specified if \code{measure="quantile"}.}
  \item{trim}{s. \code{trimmed.mean}.}
}

\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\seealso{\link{RepeatRanking}, \link{AggregateSVD}, \link{AggregatePenalty}, \link{AggregateMC}}
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
loo <- GenerateFoldMatrix(y = yy, k=1)
### Get all rankings
loor_ordT <- RepeatRanking(ordT, loo)
### aggregate rankings
agg_simple_ordT <- AggregateSimple(loor_ordT, measure ="mean")
toplist(agg_simple_ordT)
}
