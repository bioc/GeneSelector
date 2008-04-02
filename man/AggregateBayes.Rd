\name{AggregateBayes}
\alias{AggregateBayes}
\title{Bayesian aggregation of repeated rankings}
\description{
  The aggregated rank results from a posterior characteristic
  (argument \code{posteriorfun} below). The discrete
  prior is symmetrically centered around the rank obtained
  from the original dataset. The Likelihood is based on
  a normal distribution with variance \code{sigma} (s. below).
}
\usage{
AggregateBayes(RR, S, tau, sigma = c("MAD", "sd"), 
                posteriorfun = c("mode", "mean", "median", "quantile"), 
                q = NULL)
}
\arguments{
  \item{RR}{An object of class \code{RepeatRanking}.}
  \item{S}{Either an object of class \code{StabilityLm} or \code{StabilityOverlap}.}
  \item{tau}{The prior variance. Controls the confidence 
             in the rank obtained from the original dataset.\cr
             Should not be too large (<=1) in order to save
             computing time.}
  \item{sigma}{How the standard deviation for the Likelihood is to be estimated
               from the data (=ranks from perturbed datasets). \code{"MAD"}
               is a (weighted) MAD,  \code{"sd"} a (weighted) standard deviation.}
  \item{posteriorfun}{Which statistic should be applied to the posterior
                      distribution as a summary. If \code{"quantile"} is
                      chosen, then it should be specified via
                      the argument \code{q}. }
  \item{q}{The posterior quantile used as summary statistic.\cr
           Only used if \code{posteriorfun} is \code{"quantile"}}
}
\details{
 The prior has support only in the range \code{[r0-2*tau;r0+2*tau]},
 where \code{r0} is the prior mode (rank from the original dataset).\cr
 The weights for the estimation of \code{sigma} decrease linearly with
 decreasing similarity of perturbed dataset and original dataset
 as measured by Stability Measures (object \code{S}). 
}
\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, \link{GetStabilityOverlap},
         \link{AggregateSimple}}
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
agg_ordT <- AggregateBayes(loor_ordT, stab_overlap, tau=1)
}

