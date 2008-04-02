### filename: AggregateSimple.r
### Title: Aggregation using simple statistics.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 31.8.2007
### date(s) of updates: 3.9.2007
#
### Brief description:
#
#  Aggregation is based on simple (weighted) descriptive
#  statistics; weights depend on stability measures.
#
#
### Further comments and notes:
#
#   s. also: AggregateBayes.r
#   much simpler/faster and probably also more stable.
#
###**************************************************************************###

################################################################################

setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))

setClass(Class="RepeatRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))

################################################################################


setGeneric("AggregatePenalty", function(RR, lambda = NULL, k=5, theta = 50,
            estimator = c("var", "mad", "iqr", "residuals"), ...)
            standardGeneric("AggregatePenalty"))


setMethod("AggregatePenalty", signature(RR="RepeatRanking"),
            function(RR,  lambda = NULL, k=5, theta = 50,
                     estimator = c("var", "mad", "iqr", "residuals"), ...) {
if(!is.null(lambda) && lambda < 0) stop("'lambda' must be strictly greater than 0 \n")
if(is.null(lambda)) lambda <- (k-1)/(theta - 1)
if(lambda < 0 || !is.finite(lambda)) stop("Invalid values specified for 'k' or 'theta' \n")
r0 <- RR@original@ranking
ranking <- RR@rankings
ord0 <- match(1:nrow(ranking), r0)
estimator <- match.arg(estimator)
if(is.element(estimator, c("var", "mad", "iqr")))
varhat <- GeneSelector:::variance(RR, estimator,...)
else{
  stab <- GetStabilityLm(RR, ...)
  varhat <- stab@residuals.unscaled
}
agg <- ((1:nrow(ranking))*(1+lambda*varhat))[ord0]
names(agg) <- NULL
new("AggregatedRanking", posterior=NA, summary=agg, pval=RR@original@pval,
     type="penalty", fun=estimator, method = RR@original@method)
 }
)

