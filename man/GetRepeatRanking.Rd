\name{GetRepeatRanking}
\alias{GetRepeatRanking}
\title{Repeat the ranking procedure for altered data sets}
\description{
 Altered data sets are typically prepared by calls to
 \link{GenerateFoldMatrix} or \link{GenerateBootMatrix}. The
 ranking procedure is then repeated for each of these
 new 'artificial' altered data sets. One major goal of this
 procedure is to examine the stability of the results obtained
 with the original dataset. 
}
\usage{GetRepeatRanking(R, P, scheme=c("Subsampling", "Labelexchange"), iter=10, 
                              varlist = list(genewise=FALSE, factor=1/5), ...)}

\arguments{
  \item{R}{The original ranking, represented by an object of
           class \link{GeneRanking}.}
  \item{P}{An object of class \link{FoldMatrix} or \link{BootMatrix} as generated
by \link{GenerateFoldMatrix} or \link{GenerateBootMatrix}, respectively.\cr
           Can also be \code{missing}. In this case, the original dataset
           is perturbed by adding gaussian noise, s. argument
           \code{varlist}.}
  \item{scheme}{Used only if \code{P} is a \code{Foldmatrix}.
                Can be \code{"Subsampling"} or \code{"Labelexchange"}.
                'Subsampling' means that observations are removed
                as determined by the slot \code{foldmatrix}. 
                'Labelexchange' means that those observations
                which would be removed are instead kept in the sample,
                but are assigned to the opposite class.}
  \item{iter}{Used only if \code{P} is missing, specifying the
              number of different noise-perturbed datasets to be created.
              Per default, the number of iterations is 10.}
  \item{varlist}{Used only if \code{P} is missing. A list with two components
        (\code{genewise}, a logical and \code{frac}, a positive real number),
        both controlling the variance of the added noise. If 
         \code{genewise=FALSE} (default) then the noise has the same variance for
		all genes: it is estimated by pooled variance estimation from the
		original data set. Otherwise, the variance of the noise is different
		for each gene and estimated genewise from the original data set. \code{frac}
         is the fraction of the variance of the estimated variance(s)
         to be used as the variance of the added noise. The default value is \code{1/5}
         and is usually clearly smaller than 1.}
  
  \item{\dots}{Further arguments to be passed to the Ranking method
               from which rankings are generated.}
}

\value{An object of class \link{RepeatRanking}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GeneRanking}, \link{RepeatRanking},
         \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, \link{RankingLimma}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}, \link{RankingGap}}
\keyword{univar}         
\examples{
 ## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### Get ranking for the original data set, with the ordinary t-statistic
ordT <- RankingTstat(xx, yy, type="unpaired")
### Generate the leave-one-out / exchange-one-label matrix
loo <- GenerateFoldMatrix(xx, yy, k=1)
### Repeat the ranking with the t-statistic, using the leave-one-out scheme
loor_ordT <- GetRepeatRanking(ordT, loo)
### .. or the label exchange scheme
ex1r_ordT <- GetRepeatRanking(ordT, loo, scheme = "Labelexchange")
### Generate the bootstrap matrix
boot <- GenerateBootMatrix(xx, yy, maxties=3, minclassize=5, repl=30)
### Repeat ranking with the t-statistic for bootstrap replicates
boot_ordT <- GetRepeatRanking(ordT, boot)
### Repeat the ranking procedure for an altered data set with added noise
noise_ordT <- GetRepeatRanking(ordT, varlist=list(genewise=TRUE, factor=1/10))
}
