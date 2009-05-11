\name{AggregateSVD}
\alias{AggregateSVD}
\alias{AggregateSVD-methods}
\alias{AggregateSVD,RepeatedRanking-method}
\title{Agregation of repeated rankings using the singular value decomposition (SVD)}
\description{
  A matrix storing all rankings is centered rowwise (=genewise), and then approximated
  using only the first singular value and the first singular vectors (s. Golub and Van Loan (1983) for details about the SVD).
  The rowwise mean vector is added afterwards, and the rowwise mean are finally used as aggregation. A weighting scheme
  giving more weight to top genes is incorporated by an (iteratively) weighted SVD, which is re-computed until convergence.
  Note that the SVD is closely related to principal component analysis, a standard tool for dimension reduction in high-dimensional
  datasets.
}
\usage{
  AggregateSVD(RR, weightscheme = c("original", "iterative"), decay = c("linear", "quadratic", "exponential"), alpha=1)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}}
  \item{weightscheme}{If \code{weightscheme = "original"}, the weights are computed according to the reference ranking \code{RR@ranking@original}.
  If \code{weightscheme = "iterative"}, the weights are intially set to 1 for all genes. After the computation of the SVD and in turn the
  first aggregation, the weightes are updated according to that aggregation. This process is repeated until convergence.}
  \item{decay}{Argument controlling the weight decay of the weights
               of the summands contributing to the stability measure. If \code{decay=linear},
               then we have weight \code{1/r} for rank \code{r},
               if \code{decay=quadratic}, then the weight is \code{1/r^2}
               and if \code{decay=quadratic}, then the weight is \code{exp(-alpha*r)}
               where \code{alpha} is a tuning parameter, specified via
               the argument \code{alpha}.}
  \item{alpha}{s. \code{decay}.}
}

\value{An object of class \link{AggregatedRanking}.}
\references{Golub, G., Van Loan, C. (1983) \cr Matrix Computations \emph{John Hopkins University Press}}
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\seealso{\link{RepeatRanking}, \link{AggregateSimple}, \link{AggregatePenalty}, \link{AggregateMC}}
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
agg_svd_ordT <- AggregateSVD(loor_ordT, weightscheme = "iterative", decay = "linear")
toplist(agg_svd_ordT)
}
