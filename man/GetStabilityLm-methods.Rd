\name{GetStabilityLm-methods}
\docType{methods}
\alias{GetStabilityLm-methods}
\alias{GetStabilityLm,RepeatRanking-method}
\title{Stability measures for gene rankings}
\description{
  Assesses stability of gene rankings by regressing
  the rankings of perturbed datasets on the ranking
  of the original datasets in a weighted manner.
  The idea is that if stability is high, the resulting
  regression models fit well.
}
\section{Methods}{
The input is an object of class \code{RepeatRanking}.
\describe{

\item{RR = "RepeatRanking"}{signature 1}
}
For further argument and output information, consult \link{GetStabilityLm}.
}
\keyword{univar}

