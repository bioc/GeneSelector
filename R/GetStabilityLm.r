### filename: GetStabilityLm.r
### Title: Calculate stability measures using linear
###        regression with weights.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 23.8.2007
### date(s) of updates: 10.9.2007; 25.9.2007
#
### Brief description:
#
#  Stability measure for repeated rankings based on the linear
#  model with weights.

#
### Further comments and notes:
#
#  Preliminary version, s.GetRepeatRanking
#  Now two types of multivariate R2 can be computed.
#  The default ('wilks') is the original (more stable and faster)
#  version. The second  ('direct') is the direct generalization of R2
#  to the multivariate case. Its computation can fail if there exist
#  two response vectors that have exactly the same values. 
#  Additionally, matrix inversion is required.
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

setGeneric("GetStabilityLm", function(RR, decay = c("linear", "quadratic", "exponential"),
                             measure = c("wilks", "direct"), scheme = c("rank", "pval"), alpha=1, ...) 
                            standardGeneric("GetStabilityLm"))

setMethod("GetStabilityLm", signature(RR="RepeatRanking"),
             function(RR, decay = c("linear", "quadratic", "exponential"),
                      measure = c("wilks", "direct"), scheme = c("rank", "pval"), 
                      alpha=1, ...){
 
 measure <- match.arg(measure)
 if(!is.element(measure, eval(formals(GetStabilityLm)$measure)))
 stop("'measure' must be one of 'wilks' or 'direct' \n")                         
 decay <- match.arg(decay)
 if(!is.element(decay, eval(formals(GetStabilityLm)$decay)))
 stop("decay must be one of 'linear', 'quadratic', 'exponential' \n")
 if(alpha < 0) warning("'alpha' set to a negative value \n") 
 scheme <- match.arg(scheme)
 if(!is.element(scheme, eval(formals(GetStabilityLm)$scheme)))
 stop("scheme must be either 'rank' or 'pval' \n")
 
 Y <- RR@rankings
 lx <- nrow(Y)
 X <- as.matrix(1:lx)
 orig <- RR@original@ranking
 #X <- match(1:lx, orig)
 Y <- apply(Y, 2, function(z) match(orig,z))
 #Y <- apply(Y, 2, function(z) match(1:lx, z))
 
 ly <- ncol(Y)
 
 if(scheme == "rank")
  W <- switch(decay, linear=(1:lx)^(-1), quadratic=(1:lx)^(-2), exponential=exp(-alpha*(1:lx)))
  pval <- RR@original@pval 
  if(scheme =="pval")
  if(is.na(pval)[1]) stop("pvalues do not exitst, but scheme is 'pval'")
  else W <- switch(decay, linear=1/pval, quadratic=1/pval^2, exponential=exp(-alpha*pval))
  #W <- W/sum(W)
 
  lms  <- apply(Y, 2, function(z) lm(z ~ X, weights=W, ...))
  coefficients <- unlist(lapply(lms, function(z) coef(z)[2]))
  if(any(coefficients) < 0)
  warning("At least one coefficient is negative. \n")
  R2vec <- unlist(lapply(lms, function(z) summary(z)$r.squared))
  S <- sqrt(W)
  Yhat <- matrix(unlist(lapply(lms, predict)), nrow=lx, ncol=ly)
  Epsilon <- Y-Yhat
  ybarw <- apply(Y, 2, weighted.mean, w=W)
  #ybarw <- apply(S*Y, 2, mean)
  residuals.unscaled <- rowSums(Epsilon^2)
  Epsilon <- S*Epsilon
  residuals <- rowSums(Epsilon^2)
  Epsilon1 <- S*(Y-ybarw)
  if( measure == "wilks"){
  E <- crossprod(Epsilon)
  E1 <- crossprod(Epsilon1) 
  multivariateR2 <- 1- sum(diag(E))/sum(diag(E1))
  } 
  if ( measure == "direct"){
  uniqY <- unique(Y, MARGIN=2)
  if(ncol(uniqY) != ly)
  stop("measure 'direct' is inappropriate, because Y has duplicate columns \n") 
  termdiff <- sum(W) * ybarw %o% ybarw
  SSM <- crossprod(S*Yhat, S*Yhat) - termdiff
  SST <- crossprod(Epsilon1)
  svdcheck <- svd(SST)
  if(any(svdcheck$d < (.Machine$double.eps)^0.25))
  stop("SST matrix is computationally singular. Choose another weight scheme
        or decay parameter to avoid this problem \n")
  mR2 <- SSM %*% solve(SST)
  multivariateR2 <- sum(diag(mR2))/ly
  }
  
  #residuals <- residuals[orig] 
  #residuals.unscaled <- residuals.unscaled[orig]
  
  new("StabilityLm", coefficients=coefficients, R2vec=R2vec,
                        multivariateR2 = multivariateR2, residuals = residuals,
                        residuals.unscaled = residuals.unscaled, 
                        weightscheme = list(decay=decay, scheme=scheme, alpha=alpha))
  
  })
 

