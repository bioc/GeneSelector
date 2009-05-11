\name{HeatmapRankings}
\alias{HeatmapRankings}
\alias{HeatmapRankings-methods}
\alias{HeatmapRankings,RepeatedRanking-method}
\title{Heatmap of genes and rankings}
\description{
 Cluster genes and repeated rankings simultaneously
 based on a data matrix of ranks whose columns correspond
 to rankings and whose rows correspond to genes.
 The main goal is to compare different ranking procedures
 and to examine whether there are differences among
 them. Up to now, the Euclidean metric and complete-linkage clustering is used to generate
 the trees.}
\usage{
HeatmapRankings(RR, ind=1:100)
}

\arguments{
  \item{RR}{An object of class \link{RepeatedRanking}, usually generated from a call to \link{MergeMethods}.}
  \item{ind}{A vector of gene indices whose ranks are used
             to generate the heatmap. The number of elements
             should not be too large (not greater than 500)
             due to high time- and memory requirements.}
}

\value{A heatmap (plot).}
\references{Gentleman, R., Carey, V.J., Huber, W., Irizarry, R.A.,
            Dudoit, S. (editors), 2005.\cr
            Bioinformatics and Computational Biology Solutions
            Using R and Bioconductor, Chapter 10: Visualizing Data.
            \emph{Springer, N.Y.}}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
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
merged <- MergeMethods(list(ordinaryT, baldilongT, samT, wilc, wilcebam))
### plot the heatmap
HeatmapRankings(merged, ind=1:100)
}
