\name{HeatmapMethods}
\alias{HeatmapMethods}
\title{Heatmap of genes and ranking procedures}
\description{
 Cluster genes and ranking procedures simultanesously 
 based on a data matrix of ranks whose columns correspond
 to ranking procedures and whose rows correspond to genes.
 The main goal is to compare different ranking procedures
 and to examine whether there are big differences among 
 them. Up to now, the (totally unweighted) euclidean
 metric and complete-linkage clustering is used to generate
 the trees. It should be mentionned that this method only
 fulfills an exploratory task.  
}
\usage{
HeatmapMethods(Rlist, ind = 1:100)
}

\arguments{
  \item{Rlist}{A list of objects of class \link{GeneRanking} 
               or \link{AggregatedRanking}.}
  \item{ind}{A vector of gene indices whose ranks are used
             to generate the heatmap. The number of elements
             should not be too large (not greater than 500)
             due high time and memory requirements.}
}

\value{A heatmap (plot).}
\references{Gentleman, R., Carey, V.J, Huber, W., Irizarry, R.A, 
            Dudoit, S. (editors), 2005.\cr
            Bioinformatics and Computational Biology Solutions
            Using R and Bioconductor chapter 10, Visualizing Data.
            \emph{Springer, N.Y.}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
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
### plot the heatmap
HeatmapMethods(LL, ind=1:100)
}
