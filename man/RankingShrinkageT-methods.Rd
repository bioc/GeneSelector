\name{RankingShrinkageT-methods}
\docType{methods}
\alias{RankingShrinkageT-methods}
\alias{RankingShrinkageT,matrix,numeric-method}
\alias{RankingShrinkageT,matrix,factor-method}
\alias{RankingShrinkageT,ExpressionSet,character-method}
\title{Ranking based on the 'shrinkage t' statistic}
\description{ The shrinkage t statistic stabilizes the estimated variances
  appearing in the denominator of the statistic via a 
  James-Stein-Shrinkage approach. In this implementation,
  the shrinkage target is the median of the variances.}
\section{Methods}{
 The input (gene expression and class labels) 
 can be given in three different ways:
\describe{
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingShrinkageT}.
}
\keyword{univar}


