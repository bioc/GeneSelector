\name{PCAMethods}
\alias{PCAMethods}
\title{Principal Components Analysis for different ranking procedures}
\description{
 Aggregation over different perturbed datasets can be performed using
 one of the methods beginning with \code{Aggregate}.\cr
 For aggregation with respect to different ranking procedures, principal
 components analysis is used in the following manner: For each gene,
 the rankings obtained from, say \code{K}, procedures are stored in
 a \code{p x K} matrix, where \code{p} is the total number of genes.
 The different ranking procedures are intepreted as variables.\cr
 The first principal component is used to form an aggregated ranking.
}
\usage{
PCAMethods(Rlist)
}

\arguments{
  \item{Rlist}{A list of objects of class \link{GeneRanking}.}

}

\value{An object of class \link{PCAMethodsResult}.}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\keyword{univar}

