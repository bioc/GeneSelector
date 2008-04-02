\name{AggregatedRanking-class}
\docType{class}
\alias{AggregatedRanking-class}
\alias{AggregatedRanking}
%\alias{plot,AggregatedRanking,missing-method}
\alias{show,AggregatedRanking-method}

\title{"AggregatedRanking"}
\description{An object returned from the methods
             \code{AggregateBayes} and \code{AggregateSimple}}

\section{Slots}{
	 \describe{
    \item{\code{posterior}:}{A list of posterior distributions for
                             each gene, ordered according to gene indices,
                             not ranks. \code{NA} if \code{method="simple"}}
    \item{\code{summary}:}{A numeric vector of summary ranks for each gene.
                           \bold{In contrast} to \code{GeneRanking}, 
                           there is no ordering, i.e. the first entry
                           corresponds to the first gene \emph{index}
                           (row 1 of the expression matrix).}
    \item{\code{pval}:}{p-values from the original dataset (if exist), 
                        otherwise a vector of \code{NA}s.}
    \item{\code{type}:}{Describes the way of aggregation, either 
                        \code{"bayesian"}, \code{"simple"}, \code{"penalty"}
                        or \code{"pca"}.}
    \item{\code{fun}:}{The function used for aggregation, s.
                        \link{AggregateBayes} or \link{AggregateSimple}.}
    \item{\code{method}:}{The ranking method used.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{Use \code{show(AggregatedRanking-object)} for brief information.}
    \item{\link[plot,AggregatedRanking-method]{plot}}{
     Use \code{plot(GeneRanking-object, index=i)}
     to show the posterior distribution for gene \code{i}
     (if \code{type="bayesian"}).}
	 }
}

\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}