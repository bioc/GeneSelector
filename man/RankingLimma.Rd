\name{RankingLimma}
\alias{RankingLimma}
\alias{RankingLimma-methods}
\alias{RankingLimma,matrix,numeric-method}
\alias{RankingLimma,matrix,factor-method}
\alias{RankingLimma,ExpressionSet,character-method}
\title{Ranking based on the 'moderated' t statistic}
\description{
 The 'moderated' t statistic is based on a Bayesian hierarchical
 model which is estimated by an empirical Bayes approach (Smyth et al., 2003).
 The function is a wrapper to the functions \code{fitLm} and
 \code{eBayes} implemented in the \code{limma} package.
}
\usage{
RankingLimma(x, y, type = c("unpaired", "paired", "onesample"), gene.names = NULL, ...)
}

\arguments{
  \item{x}{A \code{matrix} of gene expression values with rows
           corresponding to genes and columns corresponding to observations or alternatively an object of class \code{ExpressionSet}.\cr
           If \code{type = paired}, the first half of the columns corresponds to 
           the first measurements and the second half to the second ones. 
           For instance, if there are 10 observations, each measured twice,
           stored in an expression matrix \code{expr}, 
           then \code{expr[,1]} is paired with \code{expr[,11]}, \code{expr[,2]}
           with \code{expr[,12]}, and so on.}
  \item{y}{If \code{x} is a matrix, then \code{y} may be
           a \code{numeric} vector or a factor with at most two levels.\cr
           If \code{x} is an \code{ExpressionSet}, then \code{y}
           is a character specifying the phenotype variable in
           the output from \code{pData}.\cr
           If \code{type = paired}, take care that the coding is
           analogously to the requirement concerning \code{x}}
  \item{type}{\describe{
                \item{"unpaired":}{two-sample test.}
                \item{"paired":}{paired test. Take care that the coding of \code{y}
                                 is correct (s. above)}
                \item{"onesample":}{\code{y} has only one level. 
                                    Test whether the true mean is different
                                    from zero.}
                   
                                   
                }}
  \item{gene.names}{An optional vector of gene names.}
  \item{\dots}{Further arguments passed to the function \code{eBayes},
               for instance the prior probability for differential
               expression. Consult the help of the \code{limma} package
               for details}
}

\value{An object of class \link{GeneRanking}.}
\references{Smyth, G. K., Yang, Y.-H., Speed, T. P. (2003).\cr 
            Statistical issues in microarray data analysis. 
            \emph{Methods in Molecular Biology 2:24, 111-136}.}
\author{Martin Slawski \cr
        Anne-Laure Boulesteix}
\seealso{
         \link{RepeatRanking}, \link{RankingTstat}, \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, 
          \link{RankingEbam}, \link{RankingWilcEbam}, \link{RankingSam}, 
          \link{RankingBstat}, \link{RankingShrinkageT}, \link{RankingSoftthresholdT}, 
          \link{RankingPermutation}}
\keyword{univar}
\examples{
### Load toy gene expression data
data(toydata)
### class labels
yy <- toydata[1,]
### gene expression
xx <- toydata[-1,]
### run RankingLimma
limma <- RankingLimma(xx, yy, type="unpaired")
}
