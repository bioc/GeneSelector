\name{GetStabilityUnion}
\alias{GetStabilityUnion}
\alias{GetStabilityUnion-methods}
\alias{GetStabilityUnion,RepeatedRanking-method}

\title{Stability measures for gene lists}
\description{
  The similarity of \emph{multiple} ordered genelists is assessed
  by counting the size of the union ('union count') for each position in the list. The higher the union count, the less stable are the
  ordered lists. Similarly to the 'overlap score' of Yang et al. (2006), Lottaz et al. (2006), we compute
  a weighted average of the union count entitled 'union score' .
}

\usage{
GetStabilityUnion(RR, decay = c("linear", "quadratic", "exponential"),
                    alpha = 1, noinformation = 0,...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatedRanking}}
  \item{decay}{Argument controlling the weight decay of the weights
               used for the computation of the union score. If \code{decay="linear"},
               then we have weight \code{1/l} for list position \code{l},
               if \code{decay="quadratic"}, then the weight is \code{1/l^2}
               and if \code{decay="exponential"}, then the weight is \code{exp(-alpha*l)}
               where \code{alpha} is a tuning parameter, specified via
               the argument \code{alpha}. Weights are used only for the union score, not
               for the union count.}
  \item{alpha}{s. \code{decay}.}
  \item{noinformation}{If \code{noinformation} is a positive integer, union count and -score in the no-information case are approximated by randomly
  generating multiple lists and computing scores \code{noinformation} times and averaging the results. Note that this procedure can be rather slow,
  depending on the number and the length of the lists to be compared.}
  \item{\dots}{Currently unused argument.}
}
\value{An object of class \link{StabilityOverlap}}

\references{Stolovitzky, G. (2003). Gene selection in microarray data: the elephant, the blind men, and our algorithms.
\emph{Current Opinion in Structural Biology 13, 370-376}.

            Jurman, G., Merler, S., Barla, A., Paoli, S., Galea, A., Furlanello, C. (2008).\cr
Algebraic stability indicators for ranked lists in molecular profiling. \emph{Bioinformatics 24, 258-264}
            
            Yang, X., Bentink, S., Scheid, S., Spang, R. (2006) \cr
            Similarities of ordered gene lists. \emph{Journal of Bioinformatics
and Computational Biology 4, 693-708}

            Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
            
\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}

\note{\itemize{\item Union count and union score are suitably normalized to fall into the unit interval for the sake of better interpretability, with zero corresponding
            to maximal instability.
            \item Note that this function yields exactly one stability score for multiple lists, as opposed to the pairwise approach in \link{GetStabilityOverlap} and \link{GetStabilityDistance}.}}

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
stab_union_ordT <- GetStabilityUnion(loor_ordT, decay="linear")
### display the result
plot(stab_union_ordT)
}
