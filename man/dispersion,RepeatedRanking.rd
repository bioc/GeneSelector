\name{dispersion,RepeatedRanking}
\alias{dispersion,RepeatedRanking-method}
\alias{dispersion}
\title{Compute genewise dispersion measures for repeated rankings}

\description{Dispersion is computed with respect to ranks, computed genewise. Three different measures
are implemented: standard deviation (\code{sd}), median absolute deviation (\code{mad}), and interquartile
range (\code{IQR}). The function is primarily intended to serve as helper function for \link{AggregatePenalty}.}

\usage{dispersion(RR, measure = c("sd", "mad", "iqr"), scheme = c("original", "symmetric", "user"), center = NULL)}
\arguments{
   \item{RR}{An object of class \code{RepeatedRanking}.}
   \item{measure}{Specifies the dispersion measure, s. description.}
   \item{scheme}{Specifies how the location parameter is computed. If \code{scheme="original"}, then the location parameter is chosen as the reference ranking
                 (slot \code{original}). If \code{scheme = "symmetric"}, then \code{original@ranking} and \code{rankings} are pooled to compute the
                 location parameter either as the joint mean (if \code{measure} = "mean") or the joint median (if \code{measure = "median"}). Alternatively,
                 the user may provide locations by using the \code{center} argument.}
   \item{center}{Location parameters to be used.
                 Used only if \code{scheme} = "user".}
 }
 
\value{A numeric vector containing the dispersion measure for each gene.}

\seealso{\link{GeneRanking}, \link{RepeatedRanking}}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\keyword{univar}
