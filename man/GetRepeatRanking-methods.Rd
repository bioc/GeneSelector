\name{GetRepeatRanking-methods}
\docType{methods}
\alias{GetRepeatRanking-methods}
\alias{GetRepeatRanking,GeneRanking,FoldMatrix,ANY,missing,missing-method}
\alias{GetRepeatRanking,GeneRanking,BootMatrix,missing,missing,missing-method}
\alias{GetRepeatRanking,GeneRanking,missing,missing,ANY,ANY-method}
\title{Repeat rankings for altered dataset}
\description{
 Altered datasets are typically prepared by calls to
 \link{GenerateFoldMatrix} or \link{GenerateBootMatrix}. The
 process of ranking is then repeated for each of these
 new 'artificial' datasets. One major goal of this
 procedure is to examine the stability of the results
 of the unchanged original dataset. 
}
\section{Methods}{
The input (an object of class \code{GeneRanking} is obligatory) 
can be given in three different ways:
\describe{

\item{R = "GeneRanking", P = "FoldMatrix"}{signature 1}

\item{R = "GeneRanking", P = "BootMatrix"}{signature 2}

\item{R = "GeneRanking", P = "missing"}{signature}
}
For further argument and output information, consult \link{GetRepeatRanking}
}

\keyword{univar}
