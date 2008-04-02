### filename: RankingBstat.r
### Title: Gene rankings via the B-statistic of Lönnstedt & Speed.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### date(s) of updates:
#
### Brief description:
#
#  Wrapper to 'stat.bayesian'  from the sma package.
#
#
### Further comments and notes:
#
#   Only method were 'unpaired' is not an admitted 'type'.
#
###**************************************************************************###

setGeneric("RankingBstat", function(x, y, type=c("paired", "onesample"),
            gene.names=NULL, ...) standardGeneric("RankingBstat"))


### signature: x=matrix, y=numeric.

setMethod("RankingBstat", signature(x="matrix", y="numeric"),
          function(x, y, type=c("paired", "onesample"), gene.names=NULL, ...){
          require(sma, quietly=TRUE)
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingBstat)$type)))
          stop("Argument 'type' must be either 'paired' or 'onesample'. \n")
          if(type == "paired"){
          taby <- table(y)
          if(length(taby) != 2)
          stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
          if(taby[1] != taby[2] || length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
          stop("Incorrect coding for type='paired'. \n")
           xx1 <- x[,1:taby[1]]
           xx2 <- x[,-c(1:taby[1])]
           diffxx <- xx2 - xx1
           outp <- stat.bayesian(diffxx,...)
           statistic <- outp$lods
           }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           outp <- stat.bayesian(x,...)
           statistic <- outp$lods
          }

          pvals <- rep(NA, nrow(x))
          ranking <- order(statistic, decreasing=TRUE)
          if(!is.null(gene.names))
            names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="Bstat")
})

### signature: x=matrix, y=factor.

setMethod("RankingBstat", signature(x="matrix", y="factor"),
          function(x, y, type=c("paired", "onesample"), gene.names=NULL, ...)
          RankingBstat(x, y=as.numeric(y), type, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingBstat", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("paired", "onesample"), gene.names=NULL, ...)
          RankingBstat(exprs(x), y=pData(x)[,y], type, gene.names, ...)
          )
