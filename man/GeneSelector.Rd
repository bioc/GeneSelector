\name{GeneSelector}
\alias{GeneSelector}
\alias{GeneSelector-methods}
\alias{GeneSelector,list-method}
\title{Select promising candidate genes}
\description{Given \code{GeneRankings} or \code{AggregatedRankings} obtained from
            several ranking procedures, the aim is to find is a unifying output. A
            threshold equal to the maximum rank/list position which is still relevant
            for the question of interest may be provided by the user, or the threshold can adaptively be determined
           via significance analysis in multiple testing procedures. Then, all genes are checked whether
            their ranks fall below this threshold \emph{consistenly} in all
            ranking procedures used. If  this holds, then the gene is selected.\cr
            A final order of the genes is defined by the following criteria
            \item{1.}{A user-defined ranking of the used ranking procedures, i.e.
                      the user decides which statistic he or she considers most important.}
            \item{2.}{'Selection', i.e. falling below the threshold.}
            \item{3.}{The obtained ranks. The rank from the most important
                      ranking procedure is considered, then that from the second
                      most important, and so on.} 
            }
\usage{
GeneSelector(Rlist, ind = NULL, indstatistic = 1:length(Rlist), 
            threshold = c("user", "BH", "qvalue", "Bonferroni", "Holm", 
            "Hochberg", "SidakSS", "SidakSD", "BY"), 
            maxrank = NULL, maxpval = 0.05)
}
\arguments{
  \item{Rlist}{A list of objects of class \code{RepeatedRanking}
               or \code{AggregatedRanking}, all based on the same data.}
  \item{ind}{Indices of genes to be considered. Defaults to all.}
  \item{indstatistic}{An index vector defining the importance of the
                     elements of \code{Rlist}. For instance, if \code{RList}
                     consists of five elements, then 
                     \code{indstatistic=c(2,4,1,3,5)} would give most
                     importance to the second element.}
  \item{threshold}{Determination of the threshold (s. description). Can be either \code{"user"},
                   in which case the threshold is specified via \code{maxrank}, or
                    an acronym for one of the following multiple testing procedures
                    (s. help file for \code{mt.rawp2adjp} in the package \code{multtest} for detailed information and references):
                    \describe{
                     \item{"BH"}{Benjamini-Hochberg procedure.}
                     \item{"qvalue"}{The q-value of Storey and Tibshirani (2003): "Statistical significance for genomewide studies".
            \emph{PNAS of the USA, 100, 9440-9445}.}
                     \item{"Bonferroni"}{Bonferroni procedure.}
                     \item{"Holm"}{Holm procedure.}
                     \item{"SidakSS"}{Sidak single-step procedure.}
                     \item{"SidakSD"}{Sidak step-down procedure.}
                     \item{"BY"}{Benjamini-Yekutieli procedure.}
                    }
                     In the latter case, the p-values of the element of
                    \code{Rlist} attributed most importance
                     (s. \code{indstatistic}) are adjusted and the
                     number of p-values falling below \code{maxpval}
                     is used as threshold rank. If the most important
                     statistic provides no p-values, then those of
                     the second most are used (if available), and so on.}
  \item{maxrank}{A positive integer specifying a user-defined threshold.}
  \item{maxpval}{Specified if \code{threshold} is \emph{not} \code{user}.}
}

\value{An object of class \link{GeneSelectorOutput}.}

\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}


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
### Detailed information about gene 4
plot(gkpval, which=4)}


