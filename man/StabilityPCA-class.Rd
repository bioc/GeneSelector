\name{StabilityPCA-class}
\docType{class}
\alias{StabilityPCA-class}
\alias{StabilityPCA}
\alias{show,StabilityPCA-method}

\title{"StabilityPCA"}
\description{An object returned from a call to \link{GetStabilityPCA}}
\section{Slots}{
	 \describe{
    \item{\code{eigenvalues}:}{Eigenvalues of the centered cross-product matrix
                               corresponding to the rank data matrix as described
                               in \link{GetStabilityPCA}, obtained
                               from principal components analysis, ordered
                               decreasingly.}
    \item{\code{measure}:}{The ratio of the maximum eigenvalue to the sum over all
                           eigenvalues, serving as stability measure.}

  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}




