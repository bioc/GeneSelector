### filename: RankingShrinkageT.r
### Title: Gene rankings via the Shrinkage t statistic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 5.9.2007
#
### Brief description:
#
#   James Stein shrinkage t-test following Opgen-Rhein & Strimmer (2007).
#
#
#
### Further comments and notes:
#
#   Produces same results as 'st'-package in one test example.
#
###**************************************************************************###

setGeneric("RankingShrinkageT", function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL,...) standardGeneric("RankingShrinkageT"))


### signature: x=matrix, y=numeric.

setMethod("RankingShrinkageT", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"), gene.names=NULL, ...){
          ly <- length(y)
          taby <- table(y)
          mode(x) <- "numeric"
          if(ly != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingShrinkageT)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample' \n")
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           m1 <- rowMeans(x)
           df1 <- ly-1
           ss1 <- rowSums((x-m1)^2)
           v1 <- ss1/df1
           wbar1 <- v1*df1/(df1+1)
           target1 <- median(v1)
           varhat1 <- apply(cbind(m1, wbar1, x), 1, function(z){
                                                   w <- (z[-c(1,2)]-z[1])^2
                                                   (df1+1)*sum((w-z[2])^2)/(df1^3)})
           lambda1 <- min(1, sum(varhat1)/sum((v1-target1)^2))
           vstar1 <- lambda1*target1 + (1-lambda1)*v1
           statistic <- m1/sqrt(vstar1/(df1+1))
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
           v1 <- ss1/df1
           wbar1 <- v1*df1/(df1+1)
           target1 <- median(v1)
           varhat1 <- apply(cbind(m1, wbar1, diffxx), 1, function(z){
                                                   w <- (z[-c(1,2)]-z[1])^2
                                                   (df1+1)*sum((w-z[2])^2)/(df1^3)})
           lambda1 <- min(1, sum(varhat1)/sum((v1-target1)^2))
           vstar1 <- lambda1*target1 + (1-lambda1)*v1
           statistic <- m1/sqrt(vstar1/(df1+1))
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
           df1 <- sum(ind)-1
           df2 <- sum(!ind)-1
           ss1 <- rowSums((xx1-m1)^2)
           v1 <- ss1/df1
           wbar1 <- v1*df1/(df1+1)
           target1 <- median(v1)
           varhat1 <- apply(cbind(m1, wbar1, xx1), 1, function(z){
                                                   w <- (z[-c(1,2)]-z[1])^2
                                                   (df1+1)*sum((w-z[2])^2)/(df1^3)})
           lambda1 <- min(1, sum(varhat1)/sum((v1-target1)^2))
           vstar1 <- lambda1*target1 + (1-lambda1)*v1
           ss2 <- rowSums((xx2-m2)^2)
           v2 <- ss2/df2
           wbar2 <- v2*df2/(df2+1)
           target2 <- median(v2)
           varhat2 <- apply(cbind(m2, wbar2, xx2), 1, function(z){
                                                   w <- (z[-c(1,2)]-z[1])^2
                                                   (df2+1)*sum((w-z[2])^2)/(df2^3)})
           lambda2 <- min(1, sum(varhat2)/sum((v2-target2)^2))
           vstar2 <- lambda2*target2 + (1-lambda2)*v2
           statistic <- (m1-m2)/sqrt(vstar1/(df1+1) + vstar2/(df2+1))
          }

          #if(pvalues) pvals <- 1 - pf(statistic^2, 1, dft)
          pvals <- rep(NA, nrow(x))
          ranking <- order(abs(statistic), decreasing=TRUE)
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
           names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="ShrinkageT")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingShrinkageT", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"), gene.names=NULL,...)
          RankingShrinkageT(x, y=as.numeric(y), type,  gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingShrinkageT", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"), gene.names=NULL,...)
          RankingShrinkageT(exprs(x), y=pData(x)[,y], type, gene.names, ...)
          )

