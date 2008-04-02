################################################################################
setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))

setClass(Class="RepeatRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))

################################################################################


setGeneric("AggregatePCA", function(RR) standardGeneric("AggregatePCA"))
            
setMethod("AggregatePCA", signature(RR="RepeatRanking"), function(RR){
  R0 <- RR@original@ranking
  rankings <- RR@rankings
  p <- nrow(rankings)
  R <- apply(rankings, 2, function(z) match(1:p, z))
  pcares <- prcomp(R, retx = TRUE, center = TRUE)
  rot <- pcares$x[,1]
  ordrot <- order(rot)
  agg <- match(1:p, ordrot)
  new("AggregatedRanking", posterior=NA, summary=agg, pval=RR@original@pval,
     type="PCA", fun="prcomp", method = RR@original@method)

 }
)
