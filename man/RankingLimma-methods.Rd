\name{RankingLimma-methods}
\docType{methods}
\alias{RankingLimma-methods}
\alias{RankingLimma,matrix,numeric-method}
\alias{RankingLimma,matrix,factor-method}
\alias{RankingLimma,ExpressionSet,character-method}
\title{Ranking based on the 'moderated' t statistic}
\description{
 The 'moderated' t statistic are based on a bayesian hierarchical
 model which is estimated by an empirical bayes approach.
 The function is a wrapper to the function \code{fitLm} and
 \code{eBayes} of the \code{limma} package.
}
\section{Methods}{
\describe{
The input (gene expression and class labels) 
can be given in three different ways:
\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingLimma}.
}
\keyword{univar}

