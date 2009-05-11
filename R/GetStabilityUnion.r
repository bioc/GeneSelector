### filename: GetStabilityUnion.r
### Title: Calculate stability measure using
###        the Union count of Jurman et al. (2008)
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 3.12.2008.
### date(s) of updates:
#
#
#
#
#
### Brief description:
#
#
#

#
### Further comments and notes:
#
#  Preliminary version, s.GetRepeatRanking
#
###**************************************************************************###

setGeneric("GetStabilityUnion", function(RR, decay = c("linear", "quadratic", "exponential"), alpha = 1, noinformation = 0,...)
                            standardGeneric("GetStabilityUnion"))

setMethod("GetStabilityUnion", signature(RR="RepeatedRanking"),
            function(RR, decay = c("linear", "quadratic", "exponential"), alpha = 1, noinformation = 0,...){

  bigR <- cbind(RR@original@ranking, RR@rankings)
  B <- ncol(bigR)
  lx <- nrow(bigR)
            
  decay <- match.arg(decay)
  if(!is.element(decay, eval(formals(GetStabilityUnion)$decay)))
  stop("'decay' must be one of 'linear', 'quadratic', 'exponential' \n")
  if(alpha < 0) warning("'alpha' set to a negative value \n")

   W <- switch(decay, linear=(1:lx)^(-1), quadratic=(1:lx)^(-2), exponential=exp(-alpha*(1:lx)))

  union <- unioncount(apply(bigR, 2, order))
  
  unionnormalized <-  1 - ((union - (1:lx))/(pmin(B*(1:lx),lx)-(1:lx)))
  
  unionscore <- cumsum(W*unionnormalized)/cumsum(W)
  if(noinformation <= 0) noinformation <- list()
  else{
    noinformation <- as.integer(noinformation)
    unionnoinfo <- matrix(nrow = lx, ncol = noinformation)
    unionscorenoinfo  <- matrix(nrow = lx, ncol = noinformation)
   for(r in 1:noinformation){
    RR <- replicate(B, sample(1:lx))
    unionr <- unioncount(RR)                    
    unionnoinfo[,r] <- 1 - ((unionr - (1:lx))/pmin(B*(1:lx),lx))
    unionscorenoinfo[,r] <- cumsum(W*unionnoinfo[r])/cumsum(W)
   }
    noinformation <- list(union = rowMeans(unionnoinfo), unionscore = rowMeans(unionscorenoinfo))
  }
 new("StabilityUnion", union = unionnormalized, unionscore = unionscore, noinformation = noinformation, decay = decay)

})

