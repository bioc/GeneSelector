setGeneric("AggregateSVD", function(RR, weightscheme = c("original", "iterative"),  decay = c("linear", "quadratic", "exponential"), alpha = 1) standardGeneric("AggregateSVD"))
            
setMethod("AggregateSVD", signature(RR="RepeatedRanking"), function(RR, weightscheme = c("original", "iterative"),
            decay = c("linear", "quadratic", "exponential"), alpha=1){

  weightscheme <- match.arg(weightscheme)
  if(!is.element(weightscheme, eval(formals(AggregateSVD)$weightscheme)))
  stop("'weightscheme' must be either 'original' or 'iterative' \n")

  decay <- match.arg(decay)
  if(!is.element(decay, eval(formals(AggregateSVD)$decay)))
  stop("decay must be one of 'linear', 'quadratic', 'exponential' \n")
  if(alpha < 0) warning("'alpha' set to a negative value \n")

  R <- cbind(RR@original@ranking, RR@rankings)
  p <- nrow(R)
  
  ### NOTE: SVD is used, centering and weighting.
  if(weightscheme == "original"){
    W <- switch(decay, linear = 1/R[,1], quadratic = 1/R[,1]^2, exponential=exp(-alpha*R[,1]))
    m <- rowMeans(R)
    Rc <- R - m
    svdR <- svd(W*Rc, nu = 1, nv = 1)
    fit <- (svdR$u  * svdR$d[1]) %*% t(svdR$v)  + m
    ranking <- rank(rowMeans(fit), ties.method = "first")
  }
  
  if(weightscheme == "iterative"){
    W <- 1
    converged <- FALSE
    m <- rowMeans(R)
    ranking <- numeric(p)
    Rc <- R - m
    while(! converged){
      svdR <- svd(W*Rc, nu = 1, nv = 1)
      fit <- (svdR$u  * svdR$d[1]) %*% t(svdR$v)  + m
      currank <- rank(rowMeans(fit), ties.method = "first")
      if(all(currank == ranking)) converged <- TRUE
      else{
             ranking <- currank
             W <- switch(decay, linear = 1/currank, quadratic = 1/currank^2, exponential=exp(-alpha*currank))
            }
    }
   }
  
  new("AggregatedRanking", ranking = ranking, type = "SVD", measure = "mean after truncated svd", method = RR@original@method)

 }
)
