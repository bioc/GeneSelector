\name{GetStabilityDistance}
\alias{GetStabilityDistance}
\alias{GetStabilityDistance-methods}
\alias{GetStabilityDistance,RepeatedRanking-method}

\title{Stability measures for gene rankings}
\description{
  The similarity of two rankings is assessed by computing a weighted distance
  measure. This function implements weighted absolute- and squared distance, a
  weighted versions of Spearman's rank correlation coefficient and Kendall's tau.
  Note that Spearman's rank correlation coefficient is not a distance measure in the
  classical sense, but can be obtained as transformation of the squared distance.
}

\usage{
GetStabilityDistance(RR, scheme = c("original", "pairwise"), measure = c("l1", "l2", "spearman", "kendall"),
decay = c("linear", "quadratic", "exponential"), alpha = 1, ...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}}
   \item{scheme}{If \code{scheme = "original"}, a reference ranking is compared with
                 alternative rankings. If \code{scheme = "pairwise"}, all possible pairs of rankings
                 are compared. The latter is normally used in the absence of a reference ranking, e.g.
                 if the agreement of different ranking procedures is of interest.}
  \item{measure}{The measure to be used. \code{measure = "l1"} computes a weighted absolute distance, \code{measure = "l2"} a weighted squared distance.
                 \code{measure = "spearman"} computes a weighted version of Spearman's rank correlation coefficient. \code{measure = "kendall"} computes
                 a weighted version of Kendall's tau. Note that, unlike in the function \code{cor} in \code{R base}, Kendall's tau ranges from \code{0} to \code{1}, and not
                 from \code{-1} to \code{1}, which is the case when \code{measure = "spearman"}. Absolute- and squared distance are suitably normalized
                 to fall into the unit interval for the sake of better interpretability, with zero corresponding to maximal instability.}
  \item{decay}{Argument controlling the weight decay of the weights
               of the summands contributing to the stability measure. If \code{decay="linear"},
               then we have weight \code{1/r} for rank \code{r},
               if \code{decay="quadratic"}, then the weight is \code{1/r^2}
               and if \code{decay="exponential"}, then the weight is \code{exp(-alpha*r)}
               where \code{alpha} is a tuning parameter, specified via
               the argument \code{alpha}.}
  \item{alpha}{s. \code{decay}.}
  \item{\dots}{Currently unused argument.}
}
\value{An object of class \link{StabilityDistance}}

\references{Jurman, G., Merler, S., Barla, A., Paoli, S., Galea, A., Furlanello, C. (2008).\cr
Algebraic stability indicators for ranked lists in molecular profiling. \emph{Bioinformatics 24, 258–264}

DeConde, R. P., Hawley, S., Falcon, S., Clegg, N., Knudsen, B., Etzioni, R.  (2006).\cr
Combining results of microarray experiments: a rank aggregation approach. \emph{Statistical Applications in Genetics and
Molecular Biology 5, 15}}
            
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}


\seealso{\link{RepeatRanking}}
\keyword{univar}
\examples{
 ## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### get ranking 
ordT <- RankingTstat(xx, yy, type="unpaired")
### Generate Leave-One-Out
loo <- GenerateFoldMatrix(y = yy, k=1)
### Repeat Ranking with t-statistic
loor_ordT <- RepeatRanking(ordT, loo)
### assess stability
stab_dis_ordT <- GetStabilityDistance(loor_ordT, scheme  = "original", measure = "spearman", decay="linear")
### for a short summary
summary(stab_dis_ordT, display = "all")
}
