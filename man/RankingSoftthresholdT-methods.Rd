\name{RankingSoftthresholdT-methods}
\docType{methods}
\alias{RankingSoftthresholdT-methods}
\alias{RankingSoftthresholdT,matrix,numeric-method}
\alias{RankingSoftthresholdT,matrix,factor-method}
\alias{RankingSoftthresholdT,ExpressionSet,character-method}
\title{Ranking via the 'soft-threshold' t-statistic}
\description{The 'soft-threshold' statistic is constructed
 using a linear regression model using the \code{L1}
 penalty (also referred to as LASSO penalty). In special
 cases (like here) the LASSO estimator can be 
 calculated analytically and is then called 'soft threshold'
 estimator.}
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}

}
For further argument and output information, consult 
\link{RankingSoftthresholdT}.
}
\keyword{univar}

