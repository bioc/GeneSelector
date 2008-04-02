################################################################################
setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))


setClass(Class="RepeatRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))

################################################################################

setGeneric("GetStabilityGLM", function(RR, decay = c("linear", "quadratic", "exponential"),
                              scheme = c("rank", "pval"), alpha=1, maxpval = 0.05,
                              method=c("raw", "BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"))
                            standardGeneric("GetStabilityGLM"))

setMethod("GetStabilityGLM", signature(RR="RepeatRanking"),
             function(RR, decay = c("linear", "quadratic", "exponential"),
                      scheme = c("rank", "pval"), alpha=1, maxpval,
                      method=c("raw", "BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BY")){
 origwarn <- options()$warn
 options(warn = -1)                     
 decay <- match.arg(decay)
 if(!is.element(decay, eval(formals(GetStabilityGLM)$decay)))
 stop("decay must be one of 'linear', 'quadratic', 'exponential' \n")
 if(alpha < 0) warning("'alpha' set to a negative value \n")
 scheme <- match.arg(scheme)
 if(!is.element(scheme, eval(formals(GetStabilityGLM)$scheme)))
 stop("scheme must be either 'rank' or 'pval' \n")

 Y <- RR@pvals
 if(is.na(Y[1,1])) stop("Ranking method used does not provide pvalues \n")
 method <- match.arg(method)
 if(!is.element(method, eval(formals(GetStabilityGLM)$method)))
 stop("Invalid 'method' specified \n")
 if(method != "raw") Yadj <- apply(Y, 2, AdjustPvalues, method = method)
 else Yadj <- Y
 Ybin <- apply(Yadj, 2, "<=", maxpval)
 mode(Ybin) <- "numeric"
 lx <- nrow(Y)
 Ycount <- rowSums(Ybin)
 ly <- ncol(Y)
 miss <- ly-Ycount
 X <- 1:lx



 if(scheme == "rank")
  W <- switch(decay, linear=X^(-1), quadratic=X^(-2), exponential=exp(-alpha*X))
  pval <- RR@original@pval
  if(scheme =="pval")
  if(is.na(pval)[1]) stop("pvalues do not exitst, but scheme is 'pval'")
  else W <- switch(decay, linear=1/pval, quadratic=1/pval^2, exponential=exp(-alpha*pval))
  #W <- W/sum(W)

  lmbin  <- apply(Ybin, 2, function(z) glm(z ~ X, weights=W, family = "binomial"))
  coefficients <- unlist(lapply(lmbin, function(z) coef(z)[2]))
  if(any(coefficients) < 0)
  warning("At least one coefficient is negative. \n")
  lmcount <- glm(cbind(Ycount, miss)~X, family = "binomial")
  devvec <- unlist(lapply(lmbin, function(z) summary(z)$deviance))
  devcount <- summary(lmcount)$deviance
  
  options(warn = origwarn)
  
  new("StabilityGLM", coefficients = coefficients, deviancevec= devvec, deviancecount = devcount,
                      weightscheme = list(decay=decay, scheme=scheme, alpha=alpha))

  })