\name{StabilityGLM-class}
\docType{class}
\alias{StabilityGLM-class}
\alias{StabilityGLM}
\alias{show,StabilityGLM-method}

\title{"StabilityGLM"}
\description{An object returned from a call to \link{GetStabilityGLM}}
\section{Slots}{
	 \describe{
    \item{\code{coefficients}:}{Slopes of the logistic regression with
                                response=1, if gene declared significant in iteration \code{b},
                                \code{b=1,...,B}, where 
                                \code{B=length(coefficients)}=no. of perturbed datasets and
                                regressor=ranks from original dataset.}
    \item{\code{deviancevec}:}{A numeric vector of the deviances belonging to the
                               regression models described shortly under \code{coefficients}.}
    \item{\code{deviancecount}:}{Deviance belonging to the logistic regression model
                                 where the responses defined under \code{coefficients}
                                 are added over the \code{B} iterations to obtain
                                 one response variable.}
    \item{\code{weightscheme}:}{A list that gives information about
                                the weighting scheme used, s. \link{GetStabilityGLM}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}}
	 }


\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}




