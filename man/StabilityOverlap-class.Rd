\name{StabilityOverlap-class}
\docType{class}
\alias{StabilityOverlap-class}
\alias{StabilityOverlap}
\alias{show,StabilityOverlap-method}
\alias{summary,StabilityOverlap-method}

\title{"StabilityOverlap"}
\description{An object returned from a call to \link{GetStabilityOverlap}.}
\section{Slots}{
	 \describe{
    \item{\code{overlap}:}{A matrix of overlap counts. The rows 
                           correspond to the position in the list (ranks),
                           columns to different perturbed datasets/replications.}
    \item{\code{scores}:}{A matrix of scores. The rows 
                           correspond to the position in the list (ranks),
                           columns to different perturbed datasets/replications.}
    \item{\code{weightscheme}:}{A list that gives information about
                                the weighting scheme used, s. \link{GetStabilityOverlap}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{summary}{Use \code{summary(object)} for summarized information, s. \link{summary,StabilityOverlap}}
    \item{plot}{Use \code{plot} for a graphical display, s. \link{plot,StabilityOverlap}}
	 }
}

\references{Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
            
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}



