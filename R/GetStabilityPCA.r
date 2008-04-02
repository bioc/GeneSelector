################################################################################

setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))


setClass(Class="RepeatRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))

################################################################################

setGeneric("GetStabilityPCA", function(RR) standardGeneric("GetStabilityPCA"))

setMethod("GetStabilityPCA", signature(RR="RepeatRanking"), function(RR){
  R0 <- RR@original@ranking
  rankings <- RR@rankings
  p <- nrow(rankings)
  R <- apply(rankings, 2, function(z) match(1:p, z))
  pcares <- prcomp(R, retx = TRUE, center = TRUE)
  eigenval <- (pcares$sdev)^2
  measure <- (eigenval/sum(eigenval))[1]
  new("StabilityPCA", eigenvalues = eigenval, measure = measure)

  })