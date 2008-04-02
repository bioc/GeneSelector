\name{toydata}
\alias{toydata}
\docType{data}
\title{Simulated gene expression dataset.}
\description{
  A matrix with rows corresponding to genes
  and colums corresponding to observations (arrays).
  The first row contains the class labels (1 and 2),
  the following 2000 rows the gene expressions.\cr
  The gene expressions were drawn from a multivariate
  normal distribution of dimension 2000 with mean
  vector zero and an unstructured simulated covariance
  matrix drawn from an Inverse Wishart distribution.\cr
  The first 40 genes are differentially expressed,
  the differences in the mean for the first class 
  were drawn from a normal distribution.
}
\usage{data(toydata)}


\examples{
data(toydata)
## extract class labels
yy <- toydata[1,]
table(yy)
## extract gene expressions
xx <- toydata[-1,]
}

\keyword{univar}

