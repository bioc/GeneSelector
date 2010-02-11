### filename: GetStabilityOverlap.r
### Title: Calculate stability measure using
###        the Overlap Score of Lottaz et al (2006).
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 24.8.2007, 27.8.2007
### date(s) of updates: 24./25.11.2008 (major)
#                       27.11.2008: argument 'scheme' introduced.
#                        3.12.2008: measure = 'union' removed, new method introduced.
#                                   now both overlap and overlap score computed
#                                   and returned automatically.
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

setGeneric("GetStabilityOverlap", function(RR, scheme = c("original", "pairwise"),
              decay = c("linear", "quadratic", "exponential"), alpha = 1, ...)
                            standardGeneric("GetStabilityOverlap"))

setMethod("GetStabilityOverlap", signature(RR="RepeatedRanking"),
            function(RR, scheme = c("original", "pairwise"),
            decay = c("linear", "quadratic", "exponential"), alpha=1, ...){
            
scheme <- match.arg(scheme)
if(!is.element(scheme, eval(formals(GetStabilityOverlap)$scheme)))
stop("'scheme' must be either 'original' or 'pairwise' \n")

#measure <- match.arg(measure)
#if(!is.element(measure, eval(formals(GetStabilityOverlap)$measure)))
#stop("'measure' must be either 'intersection' 'overlap_score' \n")

decay <- match.arg(decay)
if(!is.element(decay, eval(formals(GetStabilityOverlap)$decay)))
stop("'decay' must be one of 'linear', 'quadratic', 'exponential' \n")
if(alpha < 0) warning("'alpha' set to a negative value \n")

 lx <- length(RR@original@ranking)

 #maxscore <- 1:lx
 W <- switch(decay, linear=(1:lx)^(-1), quadratic=(1:lx)^(-2), exponential=exp(-alpha*(1:lx)))
 maxscore <- cumsum(W*(1:lx))

 if(scheme == "original"){
    R0 <- RR@original@ranking
    R <- RR@rankings
    ov <- apply(R, 2, function(z) overlap(order(R0), order(z), lx)) ### a matrix (several lists)
    overlapscore <- apply(ov, 2, function(z) cumsum(W*z))
    intersection <- 1/(1:lx) * ov
    overlapscore <- 1/maxscore * overlapscore
  }
  
 else{
  bigR <- cbind(RR@original@ranking, RR@rankings)
  Bplus <- ncol(bigR)
  pairs <- choose(Bplus, 2)
  intersection <-  overlapscore <- matrix(nrow = lx, ncol = pairs)
  k <- 1
   for(i in 1:Bplus){
     j <- i+1
    while(j <= Bplus){
    ov <- overlap(order(bigR[,i]), order(bigR[,j]), lx)
    score <- cumsum(W*ov)
    intersection[,k] <- 1/(1:lx) * ov
    overlapscore[,k] <-  1/maxscore * score
    #colnames(intersection)[k] <- colnames(overlapscore)[k] <- paste(i,j, sep = "vs.")
     j <- j+1
     k <- k+1
    }
  }
 }
 noinformation <- list(intersection = (1:lx)/lx, overlapscore =  cumsum(W*(1:lx)/lx)/cumsum(W))
                                  
 new("StabilityOverlap",  intersection = intersection, overlapscore = overlapscore,
                          noinformation = noinformation, scheme = scheme, decay = decay)

})

