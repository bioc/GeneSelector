\name{GetStabilityPCA-methods}
\docType{methods}
\alias{GetStabilityPCA-methods}
\alias{GetStabilityPCA,RepeatRanking-method}
\title{Stability measures for gene rankings}
\description{
 A principal components analysis is applied to the matrix storing
 the different rankings for each gene. The ratio of the first
 eigenvalue to the sum of all eigenvalues is used as stability measure.
 If stability is high/variability is low, then the first principal
 component will explain a large amount of the overall variation,
 leading to large first eigenvalue.
}
\section{Methods}{
The input is an object of class \code{RepeatRanking}.
\describe{
\item{RR = "RepeatRanking"}{signature 1}
}
For further argument and output information, consult \link{GetStabilityPCA}.
}
\keyword{univar}