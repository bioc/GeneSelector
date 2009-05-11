### filename: AggregatePenalty.r
### Title: Aggregation with a variance penalty.
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

setGeneric("AggregatePenalty", function(RR, dispersion = c("sd", "mad", "iqr"), center = NULL, gamma = 0.05, ...)
            standardGeneric("AggregatePenalty"))


setMethod("AggregatePenalty", signature(RR = "RepeatedRanking"),
            function(RR, dispersion = c("sd", "mad", "iqr"), center = NULL, gamma = 0.05) {

if(is.null(center)) center <- RR@original@ranking
p <- length(RR@original@ranking)                                                                                           
if(p != length(center))
stop("Length of 'center' and number of genes disagree \n")
if(any(center < 0) || any(center > p))
warning("Entries of 'center' negative or larger than the number of genes; verify correctness of 'center' \n")

dispersion <- match.arg(dispersion)
if(!is.element(dispersion, c("sd", "mad", "iqr")))
stop("Invalid choice for 'dispersion' \n")

if(gamma < 0 || gamma > 1)
stop("'gamma' may range from 0 to 1 \n")

sigmahat <- dispersion(RR = RR, measure = dispersion, scheme = "user", center = center)


new("AggregatedRanking", ranking = rank((1-gamma)*center + gamma*sigmahat), type= "penalty",
    measure = "goodness-stability score",  method = RR@original@method)
 }
)

