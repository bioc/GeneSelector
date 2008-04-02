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

setClass(Class="StabilityLm",
        representation(coefficients="numeric", R2vec="numeric",
                       multivariateR2 = "numeric", residuals = "numeric",
                       residuals.unscaled = "numeric", weightscheme = "list"))
                       

setClass(Class="StabilityOverlap",
        representation(overlap = "matrix", scores ="matrix", weightscheme = "list"))                       

################################################################################


setGeneric("AggregateSimple", function(RR, S,
            aggregatefun=c("mode", "mean", "median", "quantile"), q=NULL)
            standardGeneric("AggregateSimple"))
            

setMethod("AggregateSimple", signature(RR="RepeatRanking", S="StabilityLm"),
            function(RR, S, aggregatefun=c("mode", "mean", "median", "quantile"),
                      q=NULL){
R0 <- RR@original@ranking
R  <- RR@rankings
p <- nrow(R)
omega <- S@R2vec
fun <- match.arg(aggregatefun)
if(!is.element(fun, eval(formals(AggregateSimple)$aggregatefun)))
stop("Invalid 'aggregatefun' specified \n")
if(is.null(q) & fun  == "quantile"){
posteriorfun <- "median"
warning("aggregatefun is 'quantile', but argument 'q' is NULL; set to 0.5 \n")
}
r0 <- match(1:p, R0)
Rb <- apply(R, 2, function(z) match(1:p, z))
Rb0 <- cbind(r0, Rb)

aggregatefun <-  switch(fun,  mode=function(rb){
                                        tab <- table(rb)
                                        tabnames <- names(tab)
                                        maxtab <- tab[which.max(tab)]
                                        modes <- as.numeric(tabnames[tab == maxtab])
                                        lmodes <- length(modes)
                                        if(lmodes == 1) return(modes)
                                        else{
                                         sumw <- sapply(modes, function(z) sum(omega[rb == z]))
                                         return(modes[which.max(sumw)])
                                         }
                                        },
                                        mean=function(rb) weighted.mean(rb, c(1, omega)),
                                        median=function(rb) median(rb),
                                        quantile=function(rb) quantile(rb, q))

summary <- apply(Rb0, 1, aggregatefun)
new("AggregatedRanking", posterior=NA, summary=summary, pval=RR@original@pval, 
     type="simple", fun=fun, method = RR@original@method)
 }
)

setMethod("AggregateSimple", signature(RR="RepeatRanking", S="StabilityOverlap"),
            function(RR, S, aggregatefun=c("mode", "mean", "median", "quantile"),
                      q=NULL){
R0 <- RR@original@ranking
R  <- RR@rankings
p <-  length(R0)
omega <- S@scores[p,]
fun <- match.arg(aggregatefun)
if(!is.element(fun, eval(formals(AggregateSimple)$aggregatefun)))
stop("Invalid 'aggregatefun' specified \n")
if(is.null(q) & fun == "quantile"){
posteriorfun <- "median"
warning("aggregatefun is 'quantile', but argument 'q' is NULL; set to 0.5 \n")
}
r0 <- match(1:p, R0)
Rb <- apply(R, 2, function(z) match(1:p, z))
Rb0 <- cbind(r0, Rb)

aggregatefun <-  switch(fun,  mode=function(rb){
                                        tab <- table(rb)
                                        tabnames <- names(tab)
                                        maxtab <- tab[which.max(tab)]
                                        modes <- as.numeric(tabnames[tab == maxtab])
                                        lmodes <- length(modes)
                                        if(lmodes == 1) return(modes)
                                        else{
                                         sumw <- sapply(modes, function(z) sum(omega[rb == z]))
                                         return(modes[which.max(sumw)])
                                         }
                                        },
                                        mean=function(rb) weighted.mean(rb, c(1, omega)),
                                        median=function(rb) median(rb),
                                        quantile=function(rb) quantile(rb, q))

summary <- apply(Rb0, 1, aggregatefun)
new("AggregatedRanking", posterior=NA, summary=summary, pval=RR@original@pval, 
     type="simple", fun=fun, method = RR@original@method)
 }
)

