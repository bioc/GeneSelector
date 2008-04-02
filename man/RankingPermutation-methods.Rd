\name{RankingPermutation-methods}
\docType{methods}
\alias{RankingPermutation-methods}
\alias{RankingPermutation,matrix,numeric-method}
\alias{RankingPermutation,matrix,factor-method}
\alias{RankingPermutation,ExpressionSet,character-method}
\title{Ranking based on permutation tests.}
\description{The function is a wrapper for \code{mt.sample.teststat}
 from the package \code{multtest}. The ranking is here
 based on permutation p-values first, followed by the
 absolute value of the statistic.}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingPermutation}.
}
\keyword{univar}


