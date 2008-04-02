\name{GetStabilityOverlap}
\alias{GetStabilityOverlap}

\title{Stability measures for gene rankings}
\description{
  The similarity of two ordered genelists is assessed
  by a 'weighted cumulative sum' of the number of 
  overlaps up to a certain position in the list and is called
  'Overlap Score'. For further information, consult the reference given below.   
}

\usage{
GetStabilityOverlap(RR, decay = c("linear", "quadratic", "exponential"), 
                    scheme = c("rank", "pval"), ...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}}
   \item{scheme}{Whether ranks (\code{(scheme="rank")})
                or p-values (\code{(scheme="pval")})
                should be used as basis of weighting.}
  \item{decay}{argument controlling the weight decay for the weights 
               used in the overlap score. If \code{decay=linear},
               then we have weight \code{1/s} for rank/p-value \code{s},
               if \code{decay=quadratic}, then the weight is \code{1/s^2}
               and if \code{decay=quadratic}, then the weight is \code{exp(-s*alpha)}
               where \code{alpha} is a tuning parameter, specified via
               the argument \code{alpha}}
  \item{\dots}{Currently unused argument.}
}
\value{An object of class \link{GetStabilityOverlap}}

\references{Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
            
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\note{The computed overlap differs from the version described above in
      one point: Here, the overlap score are normalized to fall into
      the unit interval for better interpretability.}


\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}, 
          \link{RecoveryScore}, \link{GetAlpha}}
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
loo <- GenerateFoldMatrix(xx, yy, k=1)
### Repeat Ranking with t-statistic
loor_ordT <- GetRepeatRanking(ordT, loo)
### assess stability
stab_ov_ordT <- GetStabilityOverlap(loor_ordT, decay="linear")
### for a short summary
summary(stab_ov_ordT)
### for a graphical display
plot(stab_ov_ordT )
}
