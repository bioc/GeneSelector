\name{AdjustPvalues}
\alias{AdjustPvalues}
\title{P-value adjustment for mutiple testing.}
\description{
 Wrapper function for the functions \code{mt.rawp2adjp}
 from the package \code{multtest} and \code{qvalue.cal}
 from the package \code{siggenes}.  
}
\usage{
AdjustPvalues(pval, method = c("BH", "qvalue", "Bonferroni", 
               "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"))
}
\arguments{
  \item{pval}{A numeric vector of raw p-values}
  \item{method}{Several multiple testing adjustment procedures.}
}

\value{
  A numeric vector of adjusted p-values corresponding to the
  argument \code{pval}.
}
\references{Dudoit, S., Shaffer, J.P., Boldrick, J.C. (2003). \cr 
            Multiple Hypothesis Testing in Microarray Experiments
            \emph{Statistical Science, 18, 71-103}
            
            Storey, J.D., Tibshirani, R. (2003). \cr
            Statistical significance for genomewide studies.
            \emph{PNAS USA, 100, 9440-9445}
            }
\author{Martin Slawski \email{martin.slawski@campus.lmu.de} \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}}
\seealso{\link{GeneSelector}, \link{RecoveryScore}}
\keyword{univar}
\examples{
### Simulate 100 uniform random variates
rawp <- runif(100)
### Adjust with Benjamin-Hochberg procedure
adjustedp <- AdjustPvalues(rawp, method="BH")}