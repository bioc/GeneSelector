\name{GeneSelectorOutput-class}
\docType{class}
\alias{GeneSelectorOutput-class}
\alias{GeneSelectorOutput}
\alias{show,GeneSelectorOutput-method}
\alias{toplist,GeneSelectorOutput-method}
\alias{SelectedGenes-methods}
\alias{SelectedGenes}
\alias{SelectedGenes,GeneSelectorOutput-method}
\alias{plot,GeneSelectorOutput,missing-method}
\title{"GeneSelectorOutputRanking"}

\description{Object returned from a call to \code{GeneSelector}.}

\section{Slots}{
	 \describe{
    \item{\code{final}:}{Numeric vector storing the final ranks as provided by the \code{GeneSelector} procedure.}
    \item{\code{rankings}:}{Matrix of rankings used as input for the \code{GeneSelector}.}
    \item{\code{inout}:}{Matrix arranged in the same way as \code{rankmatrix},
                         but information is now binary: If \code{rankmatrix[i,j]} is smaller than the specified threshold,
                         then \code{inout[i,j]} equals \code{"+"} symbolizing selection,
                         whereas the \code{"-"} symbolizes removal.}
    \item{\code{selected}:}{The indices of those genes that fall below
                             the specified threshold. Can more conveniently be accessed
                              using \code{SelectedGenes}.}
    \item{\code{adjpval}:}{Numeric vector of adjusted p-values. \code{NA} if
                           no adjustment has been performed.}
    \item{\code{maxrank}:}{Threshold rank, either predefined by the user or
                           determined after p-value adjustment.}
    \item{\code{statistics}:}{The names of the ranking procedures used, ordered
                              according to their importance as defined by the user.}
  }
}
\section{Methods}{
  \describe{    
    \item{show}{Use \code{show(object)} for brief information.}
    \item{toplist}{Use \code{toplist(object, k=10)}  
          to display the top \code{k=10} genes according to the final ranking.}
    \item{SelectedGenes}{Use \code{SelectedGenes(object)} to show all genes that
          have been selected by the \link{GeneSelector}.}
    \item{plot}{Use \code{plot(object, which=1} to get detailed
          information about the gene with index 1, arranged in
          a pretty plot.}
   }
}

\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\keyword{univar}

