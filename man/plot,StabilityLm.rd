\name{plot,StabilityLm}
\alias{plot,StabilityLm}
\alias{plot,StabilityLm,missing-method}
\title{Visualize results from GetStabilityLm}
\description{Plots residuals from multivariate regression.
             If \code{E} is the estimated residual matrix,
             then the residual for gene \code{i} is
             \code{sum(E[i,]^2)}.}

%\usage{plot(StabilityLm-object)}
\arguments{
   \item{x}{An object of class \code{StabilityLm}}
 \item{frac}{The fraction of top genes for which the plot is done.
               Default is 1/50.}
 \item{scaled}{Should scaled residuals (according to the weights) be used ? \cr
               Default is \code{TRUE}}.
 \item{standardize}{Should residuals be transformed for unit variance and zero mean ? \cr
              Default is \code{TRUE}.}
 \item{...}{Additional arguments concerning graphical options.}
 }


\seealso{\link{GetRepeatRanking}, \link{GetStabilityLm}}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}