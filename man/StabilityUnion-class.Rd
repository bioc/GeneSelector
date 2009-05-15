\name{StabilityUnion-class}
\docType{class}
\alias{StabilityUnion-class}
\alias{StabilityUnion}
\alias{show,StabilityUnion-method}
\alias{plot,StabilityUnion,missing-method}

\title{"StabilityUnion"}
\description{An object returned from a call to \link{GetStabilityUnion}.}
\section{Slots}{
	 \describe{
    \item{\code{union}:}{A numeric vector of scores derived from union counts, normalized such that 1 is the maximum attainable value. The entries
                           correspond to the positions in the list such that the
                           \code{i}-th entry contains the score derived from the union count up to position \code{i}.}
    \item{\code{unionscore}:}{A numeric vector of union scores, arranged analogously to \code{union}.}
    \item{\code{noinformation}:}{A list containing an approximation to the expected union count and -score if \code{noinformation} in
                                 \link{GetStabilityUnion} has been set to a positive integer. Otherwise, the list is empty.}
    \item{\code{decay}}{The decay scheme used for the weights, s. \link{GetStabilityUnion}.}
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{plot}{\code{plot(object, frac=1/50, ...)} produces a graphical display
    of the union count (upper panel) and union score (lower panel) for increasing list position, ranging from \code{1} to
                \code{frac*number of genes}. The \code{...} argument may be used to modify graphical options.}
	 }
}

\references{Stolovitzky, G. (2003).\cr
Gene selection in microarray data: the elephant, the blind men and our algorithms. \emph{Current Opinion in Structural Biology 13, 370-376}

            Jurman, G., Merler, S., Barla, A., Paoli, S., Galea, A., Furlanello, C. (2008).\cr
Algebraic stability indicators for ranked lists in molecular profiling. \emph{Bioinformatics 24, 258-264}

            Yang, X., Bentink, S., Scheid, S., Spang, R. (2006) \cr
            Similarities of ordered gene lists. \emph{Journal of Bioinformatics
and Computational Biology 4, 693-708}

            Lottaz, C., Yang, X., Scheid, S., Spang, R. (2006) \cr
            OrderedList - a Bioconductor package for detecting
            similarity in ordered gene lists.
            \emph{Bioinformatics, 22, 2315-2316}}
            
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
        
\keyword{univar}



