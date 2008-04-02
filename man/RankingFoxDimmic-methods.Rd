\name{RankingFoxDimmic-methods}
\docType{methods}
\alias{RankingFoxDimmic-methods}
\alias{RankingFoxDimmic,matrix,numeric-method}
\alias{RankingFoxDimmic,matrix,factor-method}
\alias{RankingFoxDimmic,ExpressionSet,character-method}
\title{Ranking based on the t-statistic of Fox and Dimmic}
\description{Performs a two-sample bayesian t test for the following signatures:}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingFoxDimmic}.
}
\keyword{univar}
