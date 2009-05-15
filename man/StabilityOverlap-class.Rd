\name{StabilityOverlap-class}
\docType{class}
\alias{StabilityOverlap-class}
\alias{StabilityOverlap}
\alias{show,StabilityOverlap-method}
\alias{summary,StabilityOverlap-method}
\alias{plot,StabilityOverlap,missing-method}

\title{"StabilityOverlap"}
\description{An object returned from a call to \link{GetStabilityOverlap}.}
\section{Slots}{
	 \describe{
    \item{\code{intersection}:}{A matrix of intersection counts, normalized such that 1 is the maximum attainable value. The rows
                           correspond to the positions in the list such that
                           row \code{i} contains the number of common elements
                           of two lists up to position \code{i}. The columns correspond to lists,
                           either obtained from altered datasets or by applying multiple
                           ranking procedures.}
    \item{\code{overlapscore}:}{A matrix of overlap scores, arranged analogously to \code{intersection}.}
    \item{\code{noinformation}:}{A list containing two numeric vectors named \code{"intersection"} and
                                 \code{"overlap score"} containing for each list position intersection counts and overlap scores
                                 one would expect in the no-information case, corresponding to the random generation of
                                 two mutually independent lists.}
    \item{\code{scheme}}{s. \link{GetStabilityOverlap}.}
    \item{\code{decay}}{The decay scheme used for the weights, s. \link{GetStabilityOverlap}.}
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{summary}{Use \code{summary(object), measure = c("overlapscore", "intersection"), display = c("summary", "all"), position, digits = 3} for
                   condensed output. The argument \code{measure} specifies the measure one is interested in. The argument \code{display}
                   controls whether only a five-point summary is printed (\code{display = "summary"}) or whether all results
                   are printed (\code{display = "all"}). The argument \code{position} specifies the list position, e.g. if \code{position = 10}, then
                   intersection counts/overlap scores up to position 10 are summarized. The argument \code{digits} is used for number
                   formatting. Note that the output depends on \code{scheme}.}
    \item{plot}{\code{plot(object, frac=1/50, mode = c("mean", "all", "specific"), which = 1, ...)} produces a graphical display
    of the intersection count (upper panel) and the overlap score (lower panel) for increasing list position, ranging from \code{1} to
                \code{frac*number of genes}. The argument mode specifies whether this is done as average over all
                lists (\code{mode = "mean"}), for a certain alternative list \code{which} (\code{mode = "specific"}) or
                succesively in the form of separate plots for each alternative list (\code{mode = "all"}). The \code{...} argument may be used to modify
                graphical options.}
	 }
}

\references{Qiu, X., Xiao, Y., Gordon, A., Yakovlev, A., (2006) \cr
            Assessing stability of gene selection in microarray data
            analysis. \emph{BMC Bioinformatics 7, 50}

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



