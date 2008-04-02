\name{GenerateFoldMatrix-methods}
\docType{methods}
\alias{GenerateFoldMatrix-methods}
\alias{GenerateFoldMatrix,matrix,numeric-method}
\alias{GenerateFoldMatrix,matrix,factor-method}
\alias{GenerateFoldMatrix,ExpressionSet,character-method}
\title{Altered datasets via k-Jackknife or
       Label (class) exchange}
\description{
 Generates an object of class \link{FoldMatrix}
  that is then processed by \link{GetRepeatRanking}
  for the following signatures:
}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}

For further argument and output information, consult \link{GenerateFoldMatrix}.

}

\keyword{univar}