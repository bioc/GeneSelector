### filename: RankingLimma.r
### Title: Gene rankings via the moderated t-statistic of Smyth.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### date(s) of updates:
#
### Brief description:
#
#  Wrapper to 'lmFit' and 'eBayes'  from the package "limma".
#
#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("RankingLimma", function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL, ...) standardGeneric("RankingLimma"))


### signature: x=matrix, y=numeric.

setMethod("RankingLimma", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
          gene.names=NULL, ...){
          require(limma, quietly=TRUE)
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingLimma)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
           if(type == "unpaired"){
           y <- as.factor(y)
           ly <- levels(y)
           if(length(ly) != 2)
           stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
           ind  <- y==ly[1]
           ly <- levels(y)
           des <- model.matrix(~y)
           limo <- lmFit(x, des)
           outp <- eBayes(limo,...)
           statistic <- outp$t[,2]
           pvals <- outp$p.value[,2]
           }
          if(type == "paired"){
          taby <- table(y)
          if(length(taby) != 2)
          stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
          if(taby[1] != taby[2] || length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
          stop("Incorrect coding for type='paired'. \n")
           xx1 <- x[,1:taby[1]]
           xx2 <- x[,-c(1:taby[1])]
           diffxx <- xx2 - xx1
           des <- rep(1, length(y)/2)
           limo <- lmFit(diffxx, des)
           outp <- eBayes(limo,...)
           statistic <- outp$t[,1]
           pvals <- outp$p.value[,1]
           }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           des <- rep(1, length(y))
           limo <- lmFit(x, des)
           outp <- eBayes(limo,...)
           statistic <- outp$t[,1]
           pvals <- outp$p.value[,1]
          }
           
          ranking <- rank(-abs(statistic))
          if(!is.null(gene.names))
            names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="Limma")
})

### signature: x=matrix, y=factor.

setMethod("RankingLimma", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"), gene.names=NULL, ...)
          RankingLimma(x, y=as.numeric(y), type, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingLimma", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"), gene.names=NULL, ...)
          RankingLimma(exprs(x), y=pData(x)[,y], type, gene.names, ...)
          )

