\name{GetStabilityGLM-methods}
\docType{methods}
\alias{GetStabilityGLM-methods}
\alias{GetStabilityGLM,RepeatRanking-method}
\title{Stability measures for significance findings}
\description{
 Assesses the stability of the set of genes declared statistically
  significant for differential expression. To this end, p-values
  or adjusted p-values are used to generate binary response variables
  for a logistic regression model. As single covariate, the ranks
  obtained from the original dataset are used.
  Analogously to the linear model approach, weights are incorporated
  to attribute more importance to higher ranked genes.
  The deviance(s) resulting from these models are used as stability measure.
}
\section{Methods}{
The input is an object of class \code{RepeatRanking}.
\describe{

\item{RR = "RepeatRanking"}{signature 1}
}
For further argument and output information, consult \link{GetStabilityGLM}.
}
\keyword{univar}

