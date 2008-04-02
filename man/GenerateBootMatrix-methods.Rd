\name{GenerateBootMatrix-methods}
\docType{methods}
\alias{GenerateBootMatrix-methods}
\alias{GenerateBootMatrix,matrix,numeric-method}
\alias{GenerateBootMatrix,matrix,factor-method}
\alias{GenerateBootMatrix,ExpressionSet,character-method}
\title{Altered datasets via bootstrap}
\description{
 Generates an object of class \link{BootMatrix}
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

For further argument and output information, consult \link{GenerateBootMatrix}.

}

\keyword{univar}

