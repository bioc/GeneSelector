\name{MergeMethods}
\alias{MergeMethods}
\alias{MergeMethods-methods}
\alias{MergeMethods,list-method}
\title{Merge rankings obtained from different ranking procedures}
\description{Converts a list containing objects of class \link{GeneRanking} into an object of class \link{RepeatedRanking}.
 }
\usage{
MergeMethods(Rlist)
}

\arguments{
  \item{Rlist}{A list consisting of objects of class \link{GeneRanking}, obtained by
               application of different methods to the same dataset.}
}

\value{An object of class \link{RepeatedRanking}. The slot \code{original} is occupied by the first element of \code{Rlist}. Note that
all information contained in the \code{GeneRanking} objects is dropped, except for the ranks.}
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
}