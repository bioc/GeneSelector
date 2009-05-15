\name{GetStabilityOverlap}
\alias{GetStabilityOverlap}
\alias{GetStabilityOverlap-methods}
\alias{GetStabilityOverlap,RepeatedRanking-method}

\title{Stability measures for gene lists}
\description{
  The similarity of two ordered genelists is assessed
  by counting the size of the intersection ('overlap') for each position in the list and by  computing a weighted cumulative sum of the number of
  overlaps up to a position in the list ('overlap score'), as suggested by Yang et al. (2006) and Lottaz et al. (2006).
}

\usage{
GetStabilityOverlap(RR, scheme = c("original", "pairwise"), decay = c("linear", "quadratic", "exponential"),
                    alpha = 1, ...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}}
   \item{scheme}{If \code{scheme = "original"}, a reference list is compared with
                 alternative lists. If \code{scheme = "pairwise"}, all possible pairs of lists
                 are compared. The latter is normally used in the absence of a reference list, e.g.
                 if the agreement of different ranking procedures is of interest.}
  \item{decay}{Argument controlling the weight decay of the weights
               necessary for the computation of the overlap score. If \code{decay="linear"},
               then we have weight \code{1/l} for list position \code{l},
               if \code{decay="quadratic"}, then the weight is \code{1/l^2}
               and if \code{decay="exponential"}, then the weight is \code{exp(-alpha*l)}
               where \code{alpha} is a tuning parameter, specified via
               the argument \code{alpha}. Weights are used only for the overlap score, not
               for the intersection count.}
  \item{alpha}{s. \code{decay}.}
  \item{\dots}{Currently unused argument.}
}
\value{An object of class \link{StabilityOverlap}}

\references{Qiu, X., Xiao, Y., Gordon, A., Yakovlev, A. (2006) \cr
            Assessing stability of gene selection in microarray data
            analysis. \emph{BMC Bioinformatics 7, 50}
            
            Yang, X., Bentink, S., Scheid, S., Spang, R. (2006) \cr
            Similarities of ordered gene lists. \emph{Journal of Bioinformatics
and Computational Biology 4, 693-708}

            Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
            
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}

\note{\itemize{
       \item Contrary to Yang et al. (2006), Lottaz et al. (2006), we consider differential expression without distinguishing between over- and underexpression.
       \item Both intersection count and overlap score are suitably normalized to fall into the unit interval for the sake of better interpretability, with zero corresponding
            to maximal instability.
       }}


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
stab_ov_ordT <- GetStabilityOverlap(loor_ordT, scheme = "original", decay="linear")
### for a short summary
summary(stab_ov_ordT, measure = "intersection", display = "all", position = 10)
### for a graphical display
plot(stab_ov_ordT)
}
