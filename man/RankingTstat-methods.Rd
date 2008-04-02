\name{RankingTstat-methods}
\docType{methods}
\alias{RankingTstat-methods}
\alias{RankingTstat,matrix,numeric-method}
\alias{RankingTstat,matrix,factor-method}
\alias{RankingTstat,ExpressionSet,character-method}
\title{Ranking with the ordinary t statistic}
\description{
 Performs an ordinary t test for the following signatures:
}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}

For further argument and output information, consult \link{RankingTstat}.

}
\keyword{univar}

