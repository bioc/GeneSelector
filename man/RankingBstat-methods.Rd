\name{RankingBstat-methods}
\docType{methods}
\alias{RankingBstat-methods}
\alias{RankingBstat,matrix,numeric-method}
\alias{RankingBstat,matrix,factor-method}
\alias{RankingBstat,ExpressionSet,character-method}
\title{Ranking based on the B-statistic.}
\description{The B-statistic was motivated in a bayesian framework described by
 Lonnstedt and Speed. It is implemented in the package \code{sma},
 the function is just a wrapper.  }
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingBstat}.
}
\keyword{univar}


