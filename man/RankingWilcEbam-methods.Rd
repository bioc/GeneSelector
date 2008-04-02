\name{RankingWilcEbam-methods}
\docType{methods}
\alias{RankingWilcEbam-methods}
\alias{RankingWilcEbam,matrix,numeric-method}
\alias{RankingWilcEbam,matrix,factor-method}
\alias{RankingWilcEbam,ExpressionSet,character-method}
\title{Ranking based on the empirical bayes approach of Efron}
\description{
  The function is a wrapper for the function \code{wilc.ebam}
  from the package \code{siggenes} that implements an
  empirical bayes mixture model approach in combination
  with the Wilcoxon statistic.
}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingWilcEbam}.
}
\keyword{univar}
