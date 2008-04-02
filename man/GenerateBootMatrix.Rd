\name{GenerateBootMatrix}
\alias{GenerateBootMatrix}
%\alias{summary,BootMatrix-method}
\title{Altered datasets via bootstrap}
\description{
  Generates an object of class \link{BootMatrix}
  that is then processed by \link{GetRepeatRanking}
}
\usage{
GenerateBootMatrix(x, y, replicates = 50, type = c("unpaired", "paired", "onesample"), maxties = NULL, minclassize = 2, balancedclass = FALSE, balancedsample = FALSE, control)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{replicates}{Number of bootstrap replicates to be generated. Should
                    rarely exceed 50.} 
  \item{type}{One of \code{"paired", "unpaired", "onesample"}, depends
              on the type of test to be performed, s. for example
              \link{RankingTstat}.}
  \item{maxties}{The maximum number of ties allowed per observation.
                 For example, \code{maxties=2} means that no observation
                 occurs more than \code{maxties+1 = 3} times in a bootstrap
                 sample.}
  \item{minclassize}{If \code{minclassize=k} for some integer \code{k},
                     then the number of observations in each class are
                     grater then or equal to \code{minclassize} for 
                     each bootstrap sample.} 
  \item{balancedclass}{If \code{balancedclass=TRUE}, then the proportions
                       of the two classes are the same for each bootstrap
                       sample. It is a shortcut for a certain value of 
                       \code{minclasssize}. May not reasonable, if class
                       proportions are unbalanced in the original sample.}
  \item{balancedsample}{Should balanced bootstrap (s.details) be performed ?}
  \item{control}{Further control arguments concerning the generation 
                 process of the bootstrap matrix, s. \link{samplingcontrol}.}
}
\details{
  For the case that \code{balancedsample=TRUE}, all other contstraints
  as imposed by \code{maxties}, \code{minclassize} and so on are ignored.
  Balanced Bootstrap (s. reference below) means that each observation
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

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}

\seealso{\link{GenerateFoldMatrix}, \link{GetRepeatRanking}}
\keyword{univar}
\examples{
## Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### Generate Boot Matrix, maximum number of ties=3, 
### minimum classize=5, 30 replications:
boot <- GenerateBootMatrix(xx, yy, maxties=3, minclassize=5, repl=30)
}
