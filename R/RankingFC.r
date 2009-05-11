### filename: RankingFC.r
### Title: Ranking bei absolute value of log-fold change.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.8.2007
### date(s) of updates:
#
### Brief description:
#
#   Naive ranking by (log-)fold change.
#   Matrix assumed to be alread  logarithm-ed. (important !)
#

#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("RankingFC", function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL, LOG=FALSE, ...) standardGeneric("RankingFC"))

setMethod("RankingFC", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"), pvalues=TRUE,
            gene.names=NULL, LOG=FALSE){
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          if(LOG) x <- log(x) 
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingFC)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
          y <- as.factor(y)
          if(type == "unpaired"){
          if(nlevels(y) != 2)
          stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
          levy <- levels(y)[1]
          ind <- (y==levy)
          FC <- rowMeans(x[,ind])-rowMeans(x[,!ind])
          }
          if(type == "paired"){
           tab <- table(y)
           if(length(tab) != 2)
           stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
           xx1 <- x[,1:tab[1]]
           xx2 <- x[,-c(1:tab[1])]
           if(tab[1] != tab[2] || length(unique(y[1:tab[1]])) != 1 | length(unique(y[-c(1:tab[1])])) != 1)
           stop("Incorrect coding for type='paired'. \n")
           FC <- rowMeans(xx2-xx1)
          }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type hase been chosen 'onesample', but y has more than one level. \n")
           FC <- rowMeans(x)
          }
          if(pvalues) pvals <-  2*(1-pnorm(abs(FC),lower.tail=TRUE))
          else pvals <- rep(NA, nrow(x))
          statistic <- FC
          ranking <- rank(-abs(statistic))
          if(!is.null(gene.names))
            names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=y, statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="Foldchange")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingFC", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL, LOG=FALSE, ...)
          RankingFC(x, y=as.numeric(y), type, pvalues, gene.names, LOG, ...))

### signature: x=ExpressionSet, y=character.

setMethod("RankingFC", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL, LOG=FALSE, ...)
          RankingFC(exprs(x), y=pData(x)[,y], type, pvalues, gene.names, LOG, ...))
          
