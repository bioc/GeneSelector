\name{GenerateFoldMatrix}
\alias{GenerateFoldMatrix}
\title{Altered datasets via k-Jackknife or Label (class) exchange}
\description{
  Generates an object of class \link{FoldMatrix}
  that is then processed by \link{GetRepeatRanking}
}

\usage{
GenerateFoldMatrix(x, y, k = 1, replicates = ifelse(k==1, ncol(x), 10), type = c("unpaired", "paired", "onesample"), minclassize = 2, balanced = FALSE, control)
}

\arguments{
   \item{x}{A \code{matrix} of gene expression values with \emph{rows}
           corresponding to genes and columns corresponding to observations. \cr
           Can alternatively an object of class \code{ExpressionSet}.\cr
           If \code{type = paired}, the first half of the columns corresponds to 
           the first measurements and the second half to the second ones. 
           For instance, if there are 10 observations, each measured twice,
           stored in an expression matrix \code{expr}, 
           then \code{expr[,1]} is paired with \code{expr[,11]}, \code{expr[,2]}
           with \code{expr[,12]}, and so on.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifyig the phenotype variable in
           the output from \code{pData}.\cr
           If \code{type = paired}, take care that the coding is
           analogously to the requirement concerning \code{x}}.
  \item{k}{Number of observations that are removed or whose labels
           are exchanged. Label exchange means that the actual label
           is replaced by the label of the other class 
           (s. \link{GetRepeatRanking}).}
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
                       \code{minclasssize}. May not reasonable, if class
                       proportions are unbalanced.}
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

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GenerateBootMatrix}, \link{GetRepeatRanking}}
\keyword{univar}
\examples{
 ## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### Generate Leave-One-Out / Exchange-One-Label matrix
loo <- GenerateFoldMatrix(xx, yy, k=1)
### A more complex example
l3o <- GenerateFoldMatrix(xx, yy, k=3, replicates=30, minclassize=5) 
}