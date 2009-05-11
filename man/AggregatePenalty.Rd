\name{AggregatePenalty}
\alias{AggregatePenalty}
\alias{AggregatePenalty-methods}
\alias{AggregatePenalty,RepeatedRanking-method}
\title{Aggregation of repeated rankings using a variance penalty approach}
\description{
  The idea behind this form of aggregation is to find a compromise between
  quality on the one hand, represented by the list position/rank, and variability
  on the other hand. The latter is assessed by calling the function \link{dispersion}.
}
\usage{
AggregatePenalty(RR, dispersion = c("sd", "mad", "iqr"), center = NULL, gamma = 0.05,...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}.}
  \item{dispersion}{The dispersion measure to be used (s. \link{dispersion}):
                      \describe{
                      \item{"sd"}{standard deviation,}
                      \item{"mad"}{median absolute deviation,}
                      \item{"iqr"}{interquartile range.}
                      }}
  \item{center}{Optional numeric vector specifying for each gene the rank serving as center/location parameter for \code{dispersion}. If \code{center = NULL},
                the reference ranking \code{RR@original@ranking} is used.}
  \item{gamma}{As basis of the aggregated ranking, the quantity \code{(1-gamma)*center + gamma * dispersion} is used, i.e. the variability aspect
               dominates as \code{gamma} tends to one.}
  \item{...}{Further arguments passed to \link{dispersion}.}
}


\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\seealso{\link{RepeatRanking}, \link{AggregateSimple}, \link{AggregateSVD}, \link{AggregateMC}}
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
agg_pen_ordT <- AggregatePenalty(loor_ordT, dispersion = "iqr", gamma = 0.3)
toplist(agg_pen_ordT)
}


