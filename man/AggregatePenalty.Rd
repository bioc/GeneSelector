\name{AggregatePenalty}
\alias{AggregatePenalty}
\title{aggregation of repeated rankings by a variance penalty}
\description{
  The idea behind this form of aggregation is to find 'reliable'
  candidate genes, i.e. those ones that are highly ranked and little
  variable at the same time. Higher variability is stronger penalized.
}
\usage{
AggregatePenalty(RR, lambda = NULL, k=5, theta = 50,
                     estimator = c("var", "mad", "iqr", "residuals"), ...)
}

\arguments{
  \item{RR}{An object of class \code{RepeatRanking}.}
  \item{lambda}{A positive real number, quantifying the amount of variance penalty.
                Default is \code{NULL}, an alternative parametrization
                using \code{k} and \code{theta} is used.}
  \item{k}{Must be specified combined with \code{theta}, s.below. Not used
           if \code{lambda} is given.}
  \item{theta}{A pragmatic way of finding an appropriate value for \code{lambda}
               is to define some threshold rank \code{theta} that is still
               considered relevant and some \code{k >= 1} that expresses
               the impprtance of the first rank as compared to the threshold
               rank.}
  \item{estimator}{The variance estimator to be used:
                      \describe{
                      \item{"var"}{The usual variance estimator.}
                      \item{"mad"}{Squared median absolute deviation.} 
                      \item{"iqr"}{Interquartile range.}
                      \item{"residuals"}{Residuals from a multivariate regression,
                                         s. \link{StabilityLm}}
                      }}
  \item{...}{Further arguments passed to \link{variance,RepeatRanking-method}.}
}


\value{An object of class \link{AggregatedRanking}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\link{AggregateSimple}, \link{AggregateBayes}, \link{AggregatePenalty},
         \link{AggregatePCA}}
\keyword{univar}

