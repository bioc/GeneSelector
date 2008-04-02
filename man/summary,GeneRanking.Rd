\name{summary,GeneRanking}
\alias{summary,GeneRanking}

\title{Summarize Gene Rankings}

\description{
  Returns a five-point summary 
  (minimum, lower and upper quartile, mean, median and maximum) 
  of statistics and p-values from an object of class
  \code{GeneRanking} arranged a two-column table.
  The second column (p-values) can be \code{NA}
  in the case that p-values have not been
  computed.   
}

%\usage{summary(GeneRanking-object)}
\arguments{
 \item{object}{An object of class \code{GeneRanking}} 
 \item{...}{Additional arguments passed to the standard \code{summary} function.}
 }
\seealso{\link{GeneRanking}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
        
\keyword{univar}

