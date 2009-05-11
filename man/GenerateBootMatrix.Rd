\name{GenerateBootMatrix}
\alias{GenerateBootMatrix}
\alias{GenerateBootMatrix-methods}
\alias{GenerateBootMatrix,missing,numeric-method}
\alias{GenerateBootMatrix,missing,factor-method}
\alias{GenerateBootMatrix,ExpressionSet,character-method}
%\alias{summary,BootMatrix-method}
\title{Altered datasets via bootstrap}
\description{
  Generates an object of class \link{BootMatrix} to be used for \link{RepeatRanking}.
}
\usage{
GenerateBootMatrix(x, y, replicates = 50, type = c("unpaired", "paired", "onesample"), maxties = NULL, minclassize = 2, balancedclass = FALSE, balancedsample = FALSE, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Only needed if \code{y} is  stored within an \code{ExpressionSet}.}
  \item{y}{\code{y} may be a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in \code{pData}.\cr
           If \code{type = "paired"}, take care that the coding is correct.}
  \item{replicates}{Number of bootstrap replicates to be generated.}
  \item{type}{One of \code{"paired", "unpaired", "onesample"}, depends
              on the type of test to be performed, s. for example
              \link{RankingTstat}.}
  \item{maxties}{The maximum number of ties allowed per observation.
                 For example, \code{maxties=2} means that no observation
                 occurs more than \code{maxties+1 = 3} times per bootstrap
                 sample.}
  \item{minclassize}{If \code{minclassize=k} for some integer \code{k},
                     then the number of observations in each class are
                     grater then or equal to \code{minclassize} for 
                     each bootstrap sample.} 
  \item{balancedclass}{If \code{balancedclass=TRUE}, then the proportions
                       of the two classes are the same for each bootstrap
                       sample. It is a shortcut for a certain value of 
                       \code{minclasssize}. May not be reasonable if class
                       proportions are unbalanced in the original sample.}
  \item{balancedsample}{Should balanced bootstrap (s. details) be performed ?}
  \item{control}{Further control arguments concerning the generation 
                 process of the bootstrap matrix, s. \link{samplingcontrol}.}
}
\details{
  For the case that \code{balancedsample=TRUE}, all other constraints
  as imposed by \code{maxties}, \code{minclassize} and so on are ignored.
  Balanced bootstrap (s. reference below) means that each observation
  occurs equally frequently (with respect to all bootstrap replications).
}

\note{
 No bootstrap sample will occur more than once, i.e. each replication
 is unique.
}

\section{warning}{If the generation process (partially) fails, try to
         reduce the constraints or change the argument \code{control}.}

\value{An object of class \code{BootMatrix}}

\references{Davison, A.C., Hinkley, D.V. (1997) \cr 
            Bootstrap Methods and their Application.
            \emph{Cambridge University Press}}

\author{Martin Slawski \cr
        Anne-Laure Boulesteix}

\seealso{\link{GenerateFoldMatrix}, \link{RepeatRanking}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### Generate Boot Matrix, maximum number of ties=3,
### minimum classize=5, 30 replications:
boot <- GenerateBootMatrix(y = yy, maxties=3, minclassize=5, repl=30)
}
