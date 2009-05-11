\name{AggregatedRanking-class}
\docType{class}
\alias{AggregatedRanking-class}
\alias{AggregatedRanking}
%\alias{plot,AggregatedRanking,missing-method}
\alias{show,AggregatedRanking-method}
\alias{toplist,AggregatedRanking-method}

\title{"AggregatedRanking"}
\description{An object returned from one of the methods \code{AggregateSimple}, \code{AggregatePenalty}, \code{AggregateMC} or \code{AggregateSVD}.}

\section{Slots}{
	 \describe{
    \item{\code{ranking}:}{A numeric vector of ranks after aggregation.}
    \item{\code{type}:}{The type of aggregation used (a character).}
    \item{\code{measure}:}{The quantity used as basis for aggregation (a character).}
    \item{\code{method}:}{The ranking method used originally, i.e. before aggregation.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(object)} for brief information.}
    \item{toplist}{Use \code{toplist(object, k=10)}
                    to display the top \code{k=10} genes in the \emph{aggregated} ranking.}
	 }
}

\author{Martin Slawski  \cr
        Anne-Laure Boulesteix}
\keyword{univar}