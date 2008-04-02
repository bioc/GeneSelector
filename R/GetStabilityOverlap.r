### filename: GetStabilityOverlap.r
### Title: Calculate stability measure using
###        the Overlap Score of Lottaz et al (2006).
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 24.8.2007, 27.8.2007
### date(s) of updates:
#
### Brief description:
#
#  Stability measure for repeated rankings based on the
#  Overlap score of Lottaz.

#
### Further comments and notes:
#
#  Preliminary version, s.GetRepeatRanking
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

setGeneric("GetStabilityOverlap", function(RR, decay = c("linear", "quadratic", "exponential"),
                         scheme = c("rank", "pval"),...)
                            standardGeneric("GetStabilityOverlap"))

setMethod("GetStabilityOverlap", signature(RR="RepeatRanking"),
            function(RR, decay = c("linear", "quadratic", "exponential"),
                         scheme = c("rank", "pval"), alpha=1, ...){

decay <- match.arg(decay)
if(!is.element(decay, eval(formals(GetStabilityOverlap)$decay)))
stop("decay must be one of 'linear', 'quadratic', 'exponential' \n")
if(alpha < 0) warning("'alpha' set to a negative value \n")
scheme <- match.arg(scheme)
if(!is.element(scheme, eval(formals(GetStabilityOverlap)$scheme)))
stop("scheme must be either 'rank' or 'pval' \n")

 R0 <- RR@original@ranking
 R <- RR@rankings
 lx <- length(R0)
 W <- switch(decay, linear=(1:lx)^(-1), quadratic=(1:lx)^(-2), exponential=exp(-alpha*(1:lx)))

 overlap <- apply(R, 2, function(z) overlap(R0, z, lx))
 scores <- apply(overlap, 2, function(z) cumsum(W*z))
 maxscore <- cumsum(W*(1:lx)) 
 scores <- 1/maxscore * scores
 
 new("StabilityOverlap", overlap = overlap, scores = scores,
                          weightscheme = list(decay=decay, scheme=scheme, alpha=alpha))

 }
)

