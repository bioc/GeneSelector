\name{RankingFC-methods}
\docType{methods}
\alias{RankingFC-methods}
\alias{RankingFC,matrix,numeric-method}
\alias{RankingFC,matrix,factor-method}
\alias{RankingFC,ExpressionSet,character-method}
\title{Ranking based on the (log) foldchange}
\description{Naive ranking that only considers difference in means
  without taking variances into account.}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}

For further argument and output information, consult 
\link{RankingFC}.
}
\keyword{univar}
