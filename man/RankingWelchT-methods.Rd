\name{RankingWelchT-methods}
\docType{methods}
\alias{RankingWelchT-methods}
\alias{RankingWelchT,matrix,numeric-method}
\alias{RankingWelchT,matrix,factor-method}
\alias{RankingWelchT,ExpressionSet,character-method}
\title{Ranking based on the Welch t statistic.}
\description{
 The Welch t statistic is a better alternative to the 
  'ordinary' t statistic in the two sample, unequal variances
  setting.
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
\link{RankingWelchT}.
}
\keyword{univar}

