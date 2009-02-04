\name{GeneSelector}
\alias{GeneSelector}
\title{Exclude genes from being candidates for differential expression}
\description{\code{GeneRankings} and \code{AggregatedRankings} from
            several statistics are unified. According to
            a user-defined or adaptively determined threshold via
            multiple testing procedures, all genes are checked whether
            they fall below this threshold \emph{consistenly} in all
            statistics used. If  this criterion is not met, then the
            gene is selected.\cr
            A final order of the genes is defined by the following criteria
             \describe{
              \item{1.}{A user-defined ranking of the used statistics, i.e.
                        the user decides which statistic is most important}
              \item{2.}{'Selection', i.e. falling below the threshold yes/no}
              \item{3.}{The obtained ranks. The rank from the most important
                        statistic is considered, then that from the second
                        most important, and so on.} 
             }
            }
\usage{
GeneSelector(Rlist, ind = NULL, indstatistic = 1:length(Rlist), 
            threshold = c("user", "BH", "qvalue", "Bonferroni", "Holm", 
            "Hochberg", "SidakSS", "SidakSD", "BY"), 
            maxrank = NULL, maxpval = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rlist}{A list of objects of class \code{RepeatedRanking}
               or \code{AggregatedRanking}, all based on the same data.}
  \item{ind}{Indices of genes to be considered. Defaults to all.}
  \item{indstatistic}{An index vector defining the importance of the
                     elements of \code{Rlist} (typically this is the importance of
                     the used statistics). For instance, if \code{RList}
                     consists of five elements, then 
                     \code{indstatistic=c(2,4,1,3,5)} would give most
                     importance to the second statistic.} 
  \item{threshold}{How the threshold is determined. Can be either \code{"user"}
                   (then the threshold is specified via \code{maxrank}) or
                    a multiple testing procedure (s. \link{AdjustPvalues}). 
                    In this case, the p-values of that element of 
                    \code{Rlist} attributed most importance
                     (s. \code{indtstatistic}) are adjusted and the
                     number of p-values falling below \code{maxpval}
                     is used as threshold rank. If the most important
                     statistic provides no p-values, then the ones of
                     the second most are used (if available), and so on.}
  \item{maxrank}{Specified if \code{threshold="user"}. A positive integer
                 that is regarded as threshold rank.}
  \item{maxpval}{Specified if \code{threshold} is \emph{not} \code{user}}
}

\value{An object of class \link{CombinedRanking}.}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}


\seealso{\link{GeneRanking}, \link{AggregatedRanking}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### Get Rankings from five different statistics
ordinaryT <- RankingTstat(xx, yy, type="unpaired")
baldilongT <- RankingBaldiLong(xx, yy, type="unpaired")
samT <- RankingSam(xx, yy, type="unpaired")
wilc <- RankingWilcoxon(xx, yy, type="unpaired")
wilcebam <- RankingWilcEbam(xx, yy, type="unpaired")
### form a list
LL <- list(ordinaryT, baldilongT, samT, wilc, wilcebam)
### order statistics (assign importance)
ordstat <-  c(3,4,2,1,5)
### start GeneSelector, threshold set to rank 50
gk50 <- GeneSelector(LL, indstatistic=ordstat, maxrank=50)
### start GeneSelector, using adaptive threshold based on p-values,
### here using the multiple testing procedure of Hochberg
gkpval <- GeneSelector(LL, indstatistic=ordstat, threshold = "BH", maxpval=0.05)
### show results
show(gkpval)
str(gkpval)
toplist(gkpval)
### which genes have been selected ?
SelectedGenes(gkpval)
### relative distance plot
plot(gkpval, top=5)
### Detailed information about gene 4
GeneInfoScreen(gkpval, which=4)}


