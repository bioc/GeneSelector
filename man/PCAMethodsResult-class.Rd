\name{PCAMethodsResult-class}
\docType{class}
\alias{PCAMethodsResult-class}
\alias{PCAMethodsResult}
\alias{show,PCAMethodsResult-method}

\title{"PCAMethodsResult"}
\description{An object returned from a call to \link{PCAMethods}}
\section{Slots}{
	 \describe{
    \item{\code{summary}:}{ numeric vector of summary ranks for each gene.
                           \bold{In contrast} to \code{GeneRanking},
                           there is no ordering, i.e. the first entry
                           corresponds to the first gene \emph{index}
                           (row 1 of the expression matrix).}
    \item{\code{eigenvalues}:}{Eigenvalues of the centered cross-product matrix
                               corresponding to the rank data matrix as described
                               in \link{PCAMethods}, obtained
                               from principal components analysis, ordered
                               decreasingly.}
    \item{\code{methods}:}{A character vector containing the names of the
                           ranking methods used.}

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




