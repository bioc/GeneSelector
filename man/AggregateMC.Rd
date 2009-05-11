\name{AggregateMC}
\alias{AggregateMC}
\alias{AggregateMC-methods}
\alias{AggregateMC,RepeatedRanking-method}
\title{Aggregation of repeated rankings using a Markov chain approach}
\description{
  All obtained rankings are aggregated on the basis of Markov chain model, in which
  each gene constitutes an element of the state space. For details, see DeConde et al. (2006).
}
\usage{
  AggregateMC(RR, maxrank, type=c("MC4", "MCT"), epsilon = 0.15)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}.}
  \item{maxrank}{Due to time- and memory requirements, the computation is limited to a reduced set of candidate genes. A gene is selected
  as candidate only if at least of one its ranks is smaller than or equal to \code{maxrank}. The remainder is assigned the rank
  \code{maxrank+1} as rank after aggregation.}
  \item{type}{Specifies the computation of the matrix of transition probabilities. If \code{type = "MC4"}, the transition probabilities
  are forced to be binary, while they may principally range from zero to one if \code{type = "MCT"}, see DeConde et al. (2006) for
  details.}
  \item{epsilon}{A second parameter concerning the computation of the transition matrix, necessary to guarantee
  ergodicity and hence existence of a unique stationary distribution of the Markov chain. The value
  \code{epsilon = 0.15}, \code{0 < epsilon < 1}, is recommended in DeConde et al. (2006).}
}

\value{An object of class \link{AggregatedRanking}.}
\references{DeConde, R. P., Hawley, S., Falcon, S., Clegg, N., Knudsen, B., Etzioni, R. (2006).\cr
Combining results of microarray experiments: a rank aggregation approach. \emph{Statistical Applications in Genetics and
Molecular Biology 5, 15}}
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\seealso{\link{RepeatRanking}, \link{AggregateSVD}, \link{AggregatePenalty}, \link{AggregateSimple}}
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
agg_MC_ordT <- AggregateMC(loor_ordT, type = "MCT", maxrank = 100)
toplist(agg_MC_ordT)
}
