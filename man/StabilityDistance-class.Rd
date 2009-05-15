\name{StabilityDistance-class}
\docType{class}
\alias{StabilityDistance-class}
\alias{StabilityDistance}
\alias{show,StabilityDistance-method}
\alias{summary,StabilityDistance-method}

\title{"StabilityDistance"}
\description{An object returned from a call to \link{GetStabilityDistance}.}
\section{Slots}{
	 \describe{
    \item{\code{scores}:}{A numeric vector of stability scores, with 0 corresponding to maximum instability and 1 to maximum
                          stability.}
    \item{\code{noinformation}:}{The stability score one would expect in the no-information case, which corresponds to the random generation
                                 of two mutually independent rankings.}
    \item{\code{scheme}}{s. \link{GetStabilityDistance}.}
    \item{\code{measure}}{The distance measure used, s. \link{GetStabilityDistance}.}
    \item{\code{decay}}{The decay scheme used for the weights, s. \link{GetStabilityDistance}.}
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{summary}{Use \code{summary(object), display = c("summary", "all"), digits = 3} for
                   condensed output. The argument \code{display}
                   controls whether only a five-point summary is printed (\code{display = "summary"}) or whether all results
                   are printed (\code{display = "all"}). The argument \code{digits} is used for number
                   formatting. Note that the output additionally depends on \code{scheme}.}

	 }
}

\references{Jurman, G., Merler, S., Barla, A., Paoli, S., Galea, A., Furlanello, C. (2008).\cr
Algebraic stability indicators for ranked lists in molecular profiling. \emph{Bioinformatics 24, 258-264}

DeConde, R. P., Hawley, S., Falcon, S., Clegg, N., Knudsen, B., Etzioni, R.  (2006).\cr
Combining results of microarray experiments: a rank aggregation approach. \emph{Statistical Applications in Genetics and
Molecular Biology 5, 15}}
            
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
        
\keyword{univar}



