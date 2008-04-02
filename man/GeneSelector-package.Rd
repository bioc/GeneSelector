\name{GeneSelector-package}
\alias{GeneSelector-package}
\docType{package}
\title{Excluding Genes from being candidates for differential 
       expressions by combining  various statistics and altered datasets.}
\description{
The name 'Geneselector' stands for the exclusion
             of genes that might be considered differentially
             expressed. 'Selected' genes are those present
             at the top of the list in various featured ranking
             methods (currently 15). In addition, the stability
             of the findings are checked by creating perturbed
             versions of the original dataset, e.g. by leaving
             samples, swapping class labels, generating bootstrap
             replicates or adding noise.}
\details{
\tabular{ll}{
Package: \tab GeneSelector\cr
Type: \tab Package\cr
Version: \tab 0.9.5\cr
Date: \tab 2008-31-1\cr
License: \tab GPL (version 2 or later) \cr
}
Most Important Steps for the workflow are:
\describe{
\item{1.}{Generate a Gene Ranking with \link{RankingTstat},
          \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, 
          \link{RankingLimma}, \link{RankingEbam}, \link{RankingWilcEbam},
            \link{RankingSam}, \link{RankingBstat}, \link{RankingShrinkageT},
            \link{RankingSoftthresholdT}, \link{RankingPermutation},
            \link{RankingGap}}
\item{2.}{Inspect the toplist using \code{toplist}.}
\item{3.}{Prepare altered datasets using \link{GenerateFoldMatrix} or
          \link{GenerateBootMatrix}}
\item{4.}{Get rankings for the altered datasets with \link{GetRepeatRanking}.}
\item{5.}{Assess stability of rankings using \link{GetStabilityLm},
          \link{GetStabilityOverlap}, \link{RecoveryScore}.}
\item{6.}{Aggregate different rankings with a bayesian approach
          with \link{AggregateBayes} or in a simple manner
          (\link{AggregateSimple}).}
\item{7.}{Inspect visually the similarity of methods using \link{HeatmapMethods}.}
\item{8.}{Combine everything into the \link{GeneSelector}.}
}
}
\author{Martin Slawski \email{martin.slawski@campus.lmu.de}, \cr
        Anne-Laure Boulesteix \url{http://www.slcmsr.net/boulesteix}

Maintainer: Martin Slawski \email{martin.slawski@campus.lmu.de}.
}
\keyword{univar}


