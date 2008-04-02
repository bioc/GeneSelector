\name{RankingWilcoxon-methods}
\docType{methods}
\alias{RankingWilcoxon-methods}
\alias{RankingWilcoxon,matrix,numeric-method}
\alias{RankingWilcoxon,matrix,factor-method}
\alias{RankingWilcoxon,ExpressionSet,character-method}
\title{Ranking based on the Wilcoxon statistic}
\description{The Wilcoxon statistic is a rank-based, 'distribution free'
              alternative. It is also closely related to the 'Area under
              the curve' (AUC) in the two sample case. The implementation
              is efficient, but still far slower than that of the t-statistic. }
\section{Methods}{
The input (gene expression and class labels) 
can be given in three different ways:
\describe{

\item{x = "matrix", y = "numeric"}{signature 1}

\item{x = "matrix", y = "factor"}{signature 2}

\item{x = "ExpressionSet", y = "character"}{signature 3}
}
For further argument and output information, consult 
\link{RankingWilcoxon}.
}
\keyword{univar}


