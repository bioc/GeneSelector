\name{RankingBaldiLong-methods}
\docType{methods}
\alias{RankingBaldiLong-methods}
\alias{RankingBaldiLong,matrix,numeric-method}
\alias{RankingBaldiLong,matrix,factor-method}
\alias{RankingBaldiLong,ExpressionSet,character-method}
\title{Ranking based on the t-statistic of Baldi and Long}
\description{Performs a bayesian t test for the following signatures:}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingBaldiLong}.
}
\keyword{univar}
