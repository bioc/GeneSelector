### filename: AggregateMC.r
### Title: Aggregation using stationary probabilities of Markov Chains.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 18.11.2008
### date(s) of updates:
#
### Brief description:
#
#  Aggregation is based on a proposal of Dwork et al.,
#  taken up by DeConde et al..
#
#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("AggregateMC", function(RR, maxrank, type=c("MC4", "MCT"), epsilon = 0.15)
            standardGeneric("AggregateMC"))


setMethod("AggregateMC", signature(RR="RepeatedRanking"), function(RR, maxrank, type=c("MC4", "MCT"), epsilon = 0.15){
R0 <- RR@original@ranking
R  <- RR@rankings
R <- cbind(R0, R)
p <- nrow(R)
M <- ncol(R)
type <- match.arg(type)
if(!is.element(type, eval(formals(AggregateMC)$type)))
stop("Invalid 'type' specified \n")
if(missing(maxrank))
stop("Argument 'maxrank' has to be specified \n")
maxrank <- as.integer(maxrank)
if(maxrank <= 1 || maxrank > p)
stop("'maxrank' must exceed 1 and must be smaller than the total number of genes \n")
U <- unique(as.numeric(apply(R, 2, function(z) match(1:maxrank, z))))
ll <- length(U)
if(ll > 500)
cat("The number of distinct genes in all top-'maxrank' lists exceeds 500; computation
         might be time-consuming or even infeasible \n")
L <- matrix(data = 0, nrow  = ll, ncol = ll)


if(type == "MCT"){
  storage.mode(R) <- "integer"
  L <- .C("mct", R = R, U = as.integer(U), lx = as.integer(p), ll = as.integer(ll), mm = as.integer(M), maxrank = as.integer(maxrank), L = L)$L
  N <- L + t(L)
  L <- L * 1/N
  diag(L) <- 0
  L[!is.finite(L)] <- 0.5
}
else{
   storage.mode(R) <- "integer" 
   L <- .C("mcfour", R = R, U = as.integer(U), lx = as.integer(p), ll = as.integer(ll), mm = as.integer(M), maxrank = as.integer(maxrank), L = L)$L
 }

  L <- L/ll
                                                                     diag(L) <- 1 - rowSums(L)                                        
  L <- (1-epsilon)*L + epsilon/ll
  one <- rep(1, ll)                                                                                                                    
  Pi <-  drop(t(one) %*% qr.solve(diag(1, ll) - L + 1))
  ranking <- numeric(p)
  ranking[-U] <- ll+1
  ranking[U] <- rank(-Pi, ties.method = "first")
  new("AggregatedRanking", ranking = ranking, type = type, measure = "stationary probabilities", method = RR@original@method)}

)



