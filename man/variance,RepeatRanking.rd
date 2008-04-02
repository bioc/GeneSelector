\name{variance,RepeatRanking}
\alias{variance,RepeatRanking-method}
\alias{variance}
\title{Compute genewise variances for ranks}

\description{One application of resampling methods is estimation of variance.
Here, variance refers to ranks, computed genewise. Three different measures
are implemented: ordinary variance, (squared)\code{mad} and the interquartile
range (\code{IQR})}

\usage{variance(RR, estimator = c("var", "mad", "iqr"), center = c("perturbed", "original"))}
\arguments{
   \item{RR}{An object of class \code{RepeatRanking}}
   \item{estimator}{Specifies the variance estimator, one of \code{var}
                    (usual variance estimator),
                    \code{mad} (squared median absolute deviation),
                    \code{iqr} (interquartile range)}
   \item{center}{Estimator for the location (mean) parameter to be used.
                 Can be either the rank from the original dataset or
                 the average rank among all perturbed datasets.}
 }
 
\value{A numeric vector containing the estimated variances corresponding
       to each gene, ordered according to the gene ranking performed 
       on the original dataset.}

\seealso{\link{GeneRanking}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}
