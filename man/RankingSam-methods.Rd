\name{RankingSam-methods}
\docType{methods}
\alias{RankingSam-methods}
\alias{RankingSam,matrix,numeric-method}
\alias{RankingSam,matrix,factor-method}
\alias{RankingSam,ExpressionSet,character-method}
\title{Ranking based on the SAM statistic}
\description{
 A wrapper function to the \code{samr} package.}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingSam}.
}
\keyword{univar}

