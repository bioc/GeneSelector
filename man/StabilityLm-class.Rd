\name{StabilityLm-class}
\docType{class}
\alias{StabilityLm-class}
\alias{StabilityLm}
\alias{show,StabilityLm-method}

\title{"StabilityLm"}
\description{An object returned from a call to \link{GetStabilityLm}}
\section{Slots}{
	 \describe{
    \item{\code{coefficients}:}{Slopes of the regression with
                                response=ranks from perturbed dataset \code{b},
                                \code{b=1,...,B}, where 
                                \code{B=length(coefficients)}=no. of perturbed datasets and
                                regressor=ranks from original dataset.}
    \item{\code{R2vec}:}{A numeric vector of the univariate coefficients
                         of determination (length is equal to that of \code{coefficients}).}
    \item{\code{multivariateR2}:}{The multivariate coefficient of determination.}
    \item{\code{residuals}:}{A numeric vector of multivariate residuals. 
                             If \code{E} is the estimated residual matrix,
                             then the residual for gene \code{i} is
                              \code{sum(E[i,]^2)}. Note that in contrast
                              to \code{residuals.unscaled}, the residual
                              matrix \code{E} has already been rescaled
                              with the square root of the weight matrix.}
    \item{\code{residuals.unscaled}:}{As \code{residuals}, but without re-scaling
                                      according to different weights.}
    \item{\code{weightscheme}:}{A list that gives information about
                                the weighting scheme used, s. \link{GetStabilityLm}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{plot}{Use \code{plot(object)} for a multivariate residual plot, 
            s. \link{plot,StabilityLm}}
	 }
}

\references{Mardia, K.V., Kent, J.T., Bibby, J.M. (1979). \cr
            Multivariate Analysis. \emph{Academic Press}}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}




