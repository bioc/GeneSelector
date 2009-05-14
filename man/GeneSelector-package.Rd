\name{GeneSelector-package}
\alias{GeneSelector-package}
\docType{package}
\title{Stability and aggregation of ranked gene lists}
\description{The term 'GeneSelector' refers to a filter selecting those genes
	     which are consistently identified as differentially
	     expressed using various statistical procedures.
             'Selected' genes are those present
             at the top of the list in various featured ranking
             methods (currently 14). In addition, the stability
             of the findings can be taken into account in the final
	     ranking by examining perturbed
             versions of the original data set, e.g. by leaving
             samples, swapping class labels, generating bootstrap
             replicates or adding noise.}
\details{
\tabular{ll}{
Package: \tab GeneSelector\cr
Type: \tab Package\cr
Version: \tab 1.5.1\cr
Date: \tab 2009-13-5\cr
License: \tab GPL (version 2 or later) \cr
}
Important steps of the workflow:
\describe{
\item{1.}{Generate a Gene Ranking with \link{RankingTstat},
          \link{RankingFC}, \link{RankingWelchT}, \link{RankingWilcoxon},
          \link{RankingBaldiLong}, \link{RankingFoxDimmic}, 
          \link{RankingLimma}, \link{RankingEbam}, \link{RankingWilcEbam},
            \link{RankingSam}, \link{RankingBstat}, \link{RankingShrinkageT},
            \link{RankingSoftthresholdT}, \link{RankingPermutation}}
\item{2.}{Inspect the toplist using \code{toplist}.}
\item{3.}{Prepare altered datasets using \link{GenerateFoldMatrix} or
          \link{GenerateBootMatrix}}
\item{4.}{Get rankings for the altered datasets with \link{RepeatRanking}.}
\item{5.}{Assess stability of rankings using \link{GetStabilityOverlap},
          \link{GetStabilityDistance}, \link{GetStabilityUnion}.}
\item{6.}{Aggregate different rankings with \link{AggregateSimple}, \link{AggregatePenalty}, \link{AggregateMC} or
          \link{AggregateSVD}.}
\item{7.}{Inspect the similarity of methods visually using \link{HeatmapRankings}.}
\item{8.}{Run the \link{GeneSelector}.}
}
}
\author{Martin Slawski \email{ms@cs.uni-sb.de}, \cr
        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}

Maintainer: Martin Slawski \email{ms@cs.uni-sb.de}.
}
\keyword{univar}


