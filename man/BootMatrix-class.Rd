\name{BootMatrix-class}
\docType{class}
\alias{BootMatrix-class}
\alias{BootMatrix}
\alias{show,BootMatrix-method}
\alias{summary,BootMatrix-method}

\title{"BootMatrix"}
\description{An object returned from \link{GenerateBootMatrix} and
             which is usually passed to \link{RepeatRanking}}
\section{Slots}{
	 \describe{
    \item{\code{bootmatrix}:}{A \code{matrix} whose number of columns equals
                              the number of replications and whose number of
                              rows equals  the number of observations. 
                              Each column contains the indices of those observations
                              that are elements of the corresponding bootstrap sample. Note that each observation
may be included several times in each column.}
    \item{\code{replicates}:}{The number of bootstrap replicates.}
    \item{\code{type}:}{One of \code{"unpaired", "paired", "onesample"}, s. \link{GeneRanking}.}
    \item{\code{maxties}:}{The maximum number of allowed ties, s. \link{GenerateBootMatrix}.}
    \item{\code{minclassize}:}{The minimum class size, s. \link{GenerateBootMatrix}}
    \item{\code{balancedclass}:}{Balanced classes, s. \link{GenerateBootMatrix}}
    \item{\code{balancedsample}:}{Balanced Bootstrap, \code{TRUE/FALSE}.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(BootMatrix)} for a brief information}
     \item{summary}{Use \code{summary(BootMatrix, repl=1:2)} to obtain
                    the frequencies of each observation for
                     replications \code{1} and \code{2}}
	 }
}
\references{Davison, A.C., Hinkley, D.V. (1997) \cr 
            Bootstrap Methods and their Application.
            \emph{Cambridge University Press}}


\seealso{\link{GenerateBootMatrix}, \link{GenerateFoldMatrix}, \link{RepeatRanking}}

\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\keyword{univar}




