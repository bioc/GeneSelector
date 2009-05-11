\name{GenerateFoldMatrix}
\alias{GenerateFoldMatrix}
\alias{GenerateFoldMatrix-methods}
\alias{GenerateFoldMatrix,missing,numeric-method}
\alias{GenerateFoldMatrix,missing,factor-method}
\alias{GenerateFoldMatrix,ExpressionSet,character-method}
\title{Altered datasets via k-Jackknife or label exchange}
\description{
  Generates an object of class \link{FoldMatrix}
  to be used for \link{RepeatRanking}.
}

\usage{
GenerateFoldMatrix(x, y, k = 1, replicates = ifelse(k==1, length(y), 10), type = c("unpaired", "paired", "onesample"), minclassize = 2, balanced = FALSE, control)
}

\arguments{
   \item{x}{Only needed if \code{y} is  stored within an \code{ExpressionSet}.}
  \item{y}{\code{y} may be a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.\cr
           If \code{type = "paired"}, take care that the coding is correct.}
  \item{k}{Number of observations that are removed or whose labels
           are exchanged. Label exchange means that the observed label
           is replaced by the label of the other class 
           (s. \link{RepeatRanking}).}
  \item{replicates}{Number of replications if \code{k>1}.}
  \item{type}{One of \code{"paired", "unpaired", "onesample"}, depends
              on the type of test to be performed, s. for example
              \link{RankingTstat}.}
  \item{minclassize}{If \code{minclassize=k} for some integer \code{k},
                     then the number of observations in each class are
                     grater then or equal to \code{minclassize} for 
                     each replication.} 
  \item{balanced}{If \code{balanced=TRUE}, then the proportions
                       of the two classes are (at least approximately) the same 
                        for each replication. It is a shortcut for a certain value of 
                       \code{minclasssize}. May not be reasonable if class
                       proportions in the given dataset are unbalanced in the original sample.}
  \item{control}{Further control arguments concerning the generation 
                 process of the fold matrix, s. \link{samplingcontrol}.}
}

\note{
 No jackknif-ed dataset will occur more than once, i.e. each replication
 is unique.
}

\section{warning}{If the generation process (partially) fails, try to
         reduce the constraints or change the argument \code{control}.}

\value{An object of class \link{FoldMatrix}.}

\references{Davison, A.C., Hinkley, D.V. (1997) \cr 
            Bootstrap Methods and their Application.
            \emph{Cambridge University Press}}

\author{Martin Slawski \cr
        Anne-Laure Boulesteix}

\seealso{\link{GenerateBootMatrix}, \link{RepeatRanking}}
\keyword{univar}
\examples{
 ## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### Generate Leave-One-Out / Exchange-One-Label matrix
loo <- GenerateFoldMatrix(y = yy, k=1)
### A more complex example
l3o <- GenerateFoldMatrix(y = yy, k=3, replicates=30, minclassize=5)
}