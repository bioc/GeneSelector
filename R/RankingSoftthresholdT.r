### filename: RankingSoftthresholdT.r
### Title: Gene rankings via the soft-thresholded t-statistic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### dates(s) of update: 7.9.2007 (hyperparameter can now optionally be omitted.)
#
### Brief description:
#
#   L_1 penalized t-statistic as proposed by Wu (2005),
#   R implementation adapted from Korbinian Strimmer ('st')
#
#
#
### Further comments and notes:
#
#  Hyperparameter tuning required; can therefore be quite slow.
#
###**************************************************************************###

setGeneric("RankingSoftthresholdT", function(x, y, type=c("unpaired", "paired", "onesample"),
            lambda = c("lowess", "cor", "user"), userlambda=NULL, gene.names=NULL,...) standardGeneric("RankingSoftthresholdT"))


### signature: x=matrix, y=numeric.

setMethod("RankingSoftthresholdT", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                   lambda = c("lowess", "cor", "user"), userlambda=NULL, gene.names=NULL, ...){
          ly <- length(y)
          taby <- table(y)
          mode(x) <- "numeric"
          if(ly != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingSoftthresholdT)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample' \n")
          method <- match.arg(lambda)
          if(!is.element(method, c("lowess","cor", "user")))
          stop("Argument 'lambda' must be one of 'lowess', 'cor' or 'user'")
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           m1 <- rowMeans(x)
           df1 <- ly-1
           ss1 <- rowSums((x-m1)^2)
           s1 <-  ss1/df1*(1/(df1+1))
           sd1 <- sqrt(s1)
           if(is.null(userlambda)) lambda <- getLambda(m1, sd1, method = method)
           else lambda <- userlambda
           num <- ifelse(m1 > lambda, m1 - sign(m1) * lambda, 0)
           den <- sqrt(s1 + lambda^2/df1)
           statistic <- num/den
          }
          if(type == "paired"){
          if(length(taby) != 2)
          stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
          if(taby[1] != taby[2] || length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
          stop("Incorrect coding for type='paired'. \n")
           xx1 <- x[,1:taby[1]]
           xx2 <- x[,-c(1:taby[1])]
           diffxx <- xx2 - xx1
           m1 <- rowMeans(diffxx)
           df1 <- ly/2-1
           ss1 <- rowSums((diffxx-m1)^2)
           s1 <-  ss1/df1*(1/(df1+1))
           sd1 <- sqrt(s1)
           if(is.null(userlambda)) lambda <- getLambda(m1, sd1, method = method)
           else lambda <- userlambda
           num <- ifelse(m1 > lambda, m1 - sign(m1) * lambda, 0)
           den <- sqrt(s1 + lambda^2/df1)
           statistic <- num/den
          }
          if(type == "unpaired"){
           y <- as.factor(y)
           ly <- levels(y)
           if(length(ly) != 2)
           stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
           ind  <- y==ly[1]
           xx1 <- x[,ind]
           xx2 <- x[,!ind]
           m1 <- rowMeans(xx1)
           m2 <- rowMeans(xx2)
           delta <- m1-m2
           df1 <- sum(ind)-1
           df2 <- sum(!ind)-1
           ss1 <- rowSums((xx1-m1)^2)
           ss2 <- rowSums((xx2-m2)^2)
           s12 <-  (ss1 + ss2)/(df1+df2)*(1/(df1+1) + 1/(df2+1))
           sd12 <- sqrt(s12)
           if(is.null(userlambda)) lambda <- getLambda(delta, sd12, method = method)
           else lambda <- userlambda
           num <- ifelse(abs(delta) > lambda, delta - sign(delta) * lambda, 0)
           den <- sqrt(s12 + lambda^2/(df1+df2))
           statistic <- num/den
          }
          #if(pvalues) pvals <- 1 - pf(statistic^2, 1, dft)
          pvals <- rep(NA, nrow(x))
          ranking <- rank(-abs(statistic))
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
          names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y = as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="SoftthresholdT")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingSoftthresholdT", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                   lambda = c("lowess", "cor", "user"), userlambda=NULL, gene.names=NULL,...)
          RankingSoftthresholdT(x, y=as.numeric(y), type, lambda, userlambda, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingSoftthresholdT", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                   lambda = c("lowess", "cor", "user"), userlambda=NULL, gene.names=NULL,...)
          RankingSoftthresholdT(exprs(x), y=pData(x)[,y], type, lambda, userlambda, gene.names, ...)
          )




