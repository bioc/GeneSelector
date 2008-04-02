\name{RankingEbam-methods}
\docType{methods}
\alias{RankingEbam-methods}
\alias{RankingEbam,matrix,numeric-method}
\alias{RankingEbam,matrix,factor-method}
\alias{RankingEbam,ExpressionSet,character-method}
\title{Ranking based on the empirical bayes approach of Efron}
\description{
  The approach of Efron and colleagues is based on a mixture model
  for subpopulations: genes that are differentially expressed and
  those that are not. The posterior probability for differential
  expression serves as statistic.
  The function described below is merely a wrapper for the
  function \code{z.ebam} from the package \code{siggenes}.
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
\link{RankingEbam}.
}
\keyword{univar}

