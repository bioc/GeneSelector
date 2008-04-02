### filename: AggregateBayes.r
### Title: Aggregation in a bayesian fashion.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 31.8.2007
### date(s) of updates: 3.9.2007
#
### Brief description:
#
#  Aggregation is based on a bayesian model.
#  Parameter tau controls the confidence in
#  the ranking of the original sample; should 
#  be chosen quite small 0-2 (time).
#
### Further comments and notes:
#
#   not well tested; speed depends massively on tau, which
#   can also vary with ranking/genes.
#   s. also: internals.r
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

setGeneric("AggregateBayes", function(RR, S, tau, sigma=c("MAD", "sd"),
            posteriorfun=c("mode", "mean", "median", "quantile"), q=NULL)
            standardGeneric("AggregateBayes"))

setMethod("AggregateBayes", signature(RR="RepeatRanking", S="StabilityLm"),
            function(RR, S, tau, sigma=c("MAD", "sd"),
                      posteriorfun=c("mode", "mean", "median", "quantile"),
                      q=NULL){
R0 <- RR@original@ranking
R  <- RR@rankings
p <- nrow(R)
omega <- S@R2vec
sigma <- match.arg(sigma)
if(!is.element(sigma, c("MAD", "sd")))
stop("argument 'sigma' must be either 'MAD' or 'sd' \n")
if(!is.element(length(tau), c(1, p)))
stop("length of tau must equal 1 or the number of ranks ! \n")
if(any(tau <= 0)) stop("tau must be strictly positive \n")
fun <- match.arg(posteriorfun)
if(!is.element(fun, eval(formals(AggregateBayes)$posteriorfun)))
stop("Invalid 'posteriorfun' specified \n")
if(is.null(q) & fun  == "quantile"){
fun <- "median"
warning("posteriorfun is 'quantile', but argument 'q' is NULL; set to 0.5 \n")
}
r0 <- match(1:p, R0)
Rb <- apply(R, 2, function(z) match(1:p, z))
if(length(tau) == 1 & sigma == "MAD"){
   Rb0 <- cbind(r0, Rb)
   postList <- apply(Rb0, 1, posterior1, omega=omega, tau=tau, p=p)
   }
if(length(tau) > 1 & sigma == "MAD"){
    Rb0 <- cbind(tau, r0, Rb)
    postList <- apply(Rb0, 1, posterior2, omega=omega, p=p)
    }
if(length(tau) == 1 & sigma == "sd"){
    Rb0 <- cbind(r0, Rb)
    postList <- apply(Rb0, 1, posterior3, omega=omega, tau=tau, p=p)
    }
if(length(tau) > 1 & sigma == "sd"){
    Rb0 <- cbind(tau, r0, Rb)
    postList <- apply(Rb0, 1, posterior4, omega=omega, p=p)
    }        

posteriorfun <-  switch(fun,   mode=function(l) as.numeric(names(which.max(l))),
                               mean=function(l) weighted.mean(as.numeric(names(l)), l),
                               median=function(l) as.numeric(names(which(cumsum(l)>0.5)[1])),
                               quantile=function(l) as.numeric(names(which(cumsum(l)>q)[1])))
                                        
summary <- unlist(lapply(postList, posteriorfun))

new("AggregatedRanking", posterior=postList, summary=summary, pval= RR@original@pval, 
     type="bayesian", fun=fun, method = RR@original@method)
 }
)


setMethod("AggregateBayes", signature(RR="RepeatRanking", S="StabilityOverlap"),
            function(RR, S, tau, sigma=c("MAD", "sd"),
                      posteriorfun=c("mode", "mean", "median", "quantile"),
                      q=NULL){
R0 <- RR@original@ranking
R  <- RR@rankings
p <- nrow(R)
omega  <- S@scores[p,]
sigma <- match.arg(sigma)
if(!is.element(sigma, c("MAD", "sd")))
stop("argument 'sigma' must be either 'MAD' or 'sd' \n")
if(!is.element(length(tau), c(1, p)))
stop("length of tau must equal 1 or the number of ranks ! \n")
if(any(tau <= 0)) stop("tau must be strictly positive \n")
fun <- match.arg(posteriorfun)
if(!is.element(fun, eval(formals(AggregateBayes)$posteriorfun)))
stop("Invalid 'posteriorfun' specified \n")
if(is.null(q) & fun  == "quantile"){
fun <- "median"
warning("posteriorfun is 'quantile', but argument 'q' is NULL; set to 0.5 \n")
}
r0 <- match(1:p, R0)
Rb <- apply(R, 2, function(z) match(1:p, z))
if(length(tau) == 1 & sigma == "MAD"){
   Rb0 <- cbind(r0, Rb)
   postList <- apply(Rb0, 1, posterior1, omega=omega, tau=tau, p=p)
   }
if(length(tau) > 1 & sigma == "MAD"){
    Rb0 <- cbind(tau, r0, Rb)
    postList <- apply(Rb0, 1, posterior2, omega=omega, p=p)
    }
if(length(tau) == 1 & sigma == "sd"){
    Rb0 <- cbind(r0, Rb)
    postList <- apply(Rb0, 1, posterior3, omega=omega, tau=tau, p=p)
    }
if(length(tau) > 1 & sigma == "sd"){
    Rb0 <- cbind(tau, r0, Rb)
    postList <- apply(Rb0, 1, posterior4, omega=omega, p=p)
    }        

posteriorfun <-  switch(fun,   mode=function(l) as.numeric(names(which.max(l))),
                               mean=function(l) weighted.mean(as.numeric(names(l)), l),
                               median=function(l) as.numeric(names(which(cumsum(l)>0.5)[1])),
                               quantile=function(l) as.numeric(names(which(cumsum(l)>q)[1])))                                        
summary <- unlist(lapply(postList, posteriorfun))
new("AggregatedRanking", posterior=postList, summary=summary, pval= RR@original@pval, 
      type="bayesian",  fun=fun, method = RR@original@method)

 }
)

