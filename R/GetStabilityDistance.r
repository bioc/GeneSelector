### filename: GetStabilityDistance.r
### Title: Calculate stability based on distance measures.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 25.11.2008
#
#
###**************************************************************************###
### NOTE 1: no partial list for 'kendall' and 'spearman'.
### NOTE 2: for canberra distance, weights are predefined. Actually, now,
### canberra distance corresponds to scheme = pairwise, decay = linear.
### NOTE 3: Return is a vector, not a matrix (as for GetStabilityOverlap)
################################################################################

setGeneric("GetStabilityDistance", function(RR,  scheme = c("original", "pairwise"), measure = c("l1", "l2", "spearman", "kendall"),
                               decay = c("linear", "quadratic", "exponential"), alpha=1, ...)
                            standardGeneric("GetStabilityDistance"))

setMethod("GetStabilityDistance", signature(RR="RepeatedRanking"),
             function(RR, scheme = c("original", "pairwise"), measure = c("l1", "l2", "spearman", "kendall"),
                         decay = c("linear", "quadratic", "exponential"), alpha=1, ...){

 scheme <- match.arg(scheme)
 if(!is.element(scheme, eval(formals(GetStabilityDistance)$scheme)))
 stop("'scheme' must be either 'original' or 'pairwise' \n")

 measure <- match.arg(measure)
 if(!is.element(measure, eval(formals(GetStabilityDistance)$measure)))
 stop("Invalid 'measure' specified \n")
 
 decay <- match.arg(decay)
 if(!is.element(decay, eval(formals(GetStabilityDistance)$decay)))
 stop("'decay' must be one of 'linear', 'quadratic', 'exponential' \n")
 if(alpha < 0) warning("'alpha' set to a negative value \n") 

 lx <- length(RR@original@ranking)
                                                  
 if(scheme == "original"){
    R0 <- RR@original@ranking
    R <- RR@rankings
    W <- switch(decay, linear = 1/R0, quadratic = 1/R0^2, exponential=exp(-alpha*R0))
    Wmaxscore <- switch(decay, linear = 1/(1:lx), quadratic = 1/((1:lx)^2), exponential = exp(-alpha*(1:lx)))
    if(is.element(measure, c("l1", "l2", "kendall")))
      maxscore <- switch(measure, l1 = sum(Wmaxscore*abs((1:lx)-(lx:1))),
                          l2 = sum(Wmaxscore*(((1:lx)-(lx:1))^2)),
                          kendall = .C("kendall", r = as.integer(1:lx),
                                                  rprime = as.integer(lx:1),
                                                  lx = as.integer(lx),
                                                  w = as.double(Wmaxscore),
                                                  res = double(1))$res)
      scores <- switch(measure, l1 = apply(R, 2, function(z) sum(W*abs(R0 - z))),
                                l2 = apply(R, 2, function(z) sum(W*(R0 - z)^2)),
                                spearman =  apply(R, 2, function(z) cor(sqrt(W)*R0, sqrt(W)*z, method = "pearson")),
                                kendall  =  apply(R, 2, function(z) .C("kendall", r = as.integer(R0),
                                                                             rprime = as.integer(z),
                                                                             lx = as.integer(lx),
                                                                             w = as.double(W),
                                                                             res = double(1))$res))
   }
  
  else{
     bigR <- cbind(RR@original@ranking, RR@rankings)
     Bplus <- ncol(bigR)
     pairs <- choose(Bplus, 2)
     scores <-  numeric(pairs)
     Wmaxscore <- switch(decay, linear = 1/(1:lx), quadratic = 1/((1:lx)^2), exponential=exp(-alpha*(1:lx)))
     if(is.element(measure, c("l1", "l2", "kendall")))
      maxscore <- switch(measure, l1 = sum(Wmaxscore*abs((1:lx)-(lx:1))),
                          l2 = sum(Wmaxscore*abs((1:lx)-(lx:1))^2),
                          kendall = .C("kendall", r = as.integer(1:lx),
                                                  rprime = as.integer(lx:1),
                                                  lx = as.integer(lx),
                                                  w = as.double(Wmaxscore),
                                                  res = double(1))$res)
     k <- 1
     for(i in 1:Bplus){
       j <- i+1
       while(j <= Bplus){                
        W <- switch(decay,  linear = 1/pmin(bigR[,i], bigR[,j]),
                            quadratic = 1/pmin(bigR[,i], bigR[,j]),
                            exponential = exp(-alpha*pmin(bigR[,i], bigR[,j])))

        scores[k] <- switch(measure, l1 = sum(W*abs(bigR[,i] - bigR[,j])),
                                  l2 = sum(W*(bigR[,i] - bigR[,j])^2),
                                  spearman =  cor(sqrt(W)*bigR[,i], sqrt(W)*bigR[,j], method = "pearson"),
                                  kendall  =  .C("kendall", r = as.integer(bigR[,i]),
                                                                             rprime = as.integer(bigR[,j]),
                                                                             lx = as.integer(lx),
                                                                             w = as.double(W),
                                                                             res = double(1))$res)
        names(scores)[k] <- paste(i,j, sep = "vs.")
        j <- j+1
        k <- k+1
      }
     }                  
     }

        noinformation <- switch(measure, l1 =  sum((2/lx^2)*(lx*lx*(lx+1)/2 -lx^3/3 - 0.5*lx^2 - lx/6) * Wmaxscore),
                                       l2 =  sum((2/lx^2)*(lx * (lx^3/3 + 0.5*lx^2 + lx/6) - 0.25*(lx+1)^4 + 0.5*(lx+1)^3 - 0.25*(lx+1)^2) * Wmaxscore),
                                       spearman = 0,
                                       kendall = 0.5 * maxscore)
                                                                                        
      if(measure != "spearman"){
        scores <- 1 - scores/maxscore
        noinformation <- 1-noinformation/maxscore
       } 
      
      new("StabilityDistance", scores =  scores, noinformation = noinformation, measure = measure, scheme = scheme, decay = decay)
       })                
  
 

