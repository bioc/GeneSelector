\name{samplingcontrol}
\alias{samplingcontrol}
\title{Control function}
\description{
  Normally, this function is not called. Only if warnings
  occur in \link{GenerateBootMatrix} or \link{GenerateFoldMatrix},
  try to increase \code{candreplicates} with respect to the default
  (three times the number of desired Boostrap/Jackknife-Iterations,
   s. argument \code{replicates} in \link{GenerateBootMatrix}/
   \link{GenerateFoldMatrix} or \code{maxiter}.)
}
\usage{
samplingcontrol(candreplicates, maxiter = 5)
}

\arguments{
  \item{candreplicates}{s. description}
  \item{maxiter}{s. description}
}

\value{A list used in \link{GenerateBootMatrix}/\link{GenerateFoldMatrix}.}
\keyword{univar}

