### filename: RankingGap.r
### Title: Gene rankings via a gaps between two groups.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### date(s) of updates:
#
### Brief description:
#
#   Ranking based on gaps between two groups so that they can
#   likely be discriminated. Idea stems from a genefilter function
#   (which is there regarded as unsupervised method).
#   If classes overlap, the statistic takes the value zero.
#   Generalization to one-class (distance from zero)
#   and paired two-class (distance of difference values from zeros)
#   possible.
#
#
### Further comments and notes:
#
#   Naive method.
#
###**************************************************************************###


### generic function

setGeneric("RankingGap", function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL,...) standardGeneric("RankingGap"))


### signature: x=matrix, y=numeric.

setMethod("RankingGap", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
          gene.names=NULL, ...){
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingGap)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
          if(type == "unpaired"){
           y <- as.factor(y)
           ly <- levels(y)
           if(length(ly) != 2)
           stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
           ind  <- y==ly[1]
           xx1 <- x[,ind]
           xx2 <- x[,!ind]
           min1 <- rowMin(xx1)
           max1 <- rowMax(xx1)
           min2 <- rowMin(xx2)
           max2 <- rowMax(xx2)
           statistic2 <- ifelse(min2 > max1, min2 - min1, 0)
           statistic1 <- ifelse(min1 > max2,  min1 - max2, 0)
           statistic <- pmax(statistic2, statistic1)
          }
          if(type == "paired"){
           tab <- table(y)
           if(length(tab) != 2)
           stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
           xx1 <- x[,1:tab[1]]
           xx2 <- x[,-c(1:tab[1])]
           if(tab[1] != tab[2] || length(unique(y[1:tab[1]])) != 1 | length(unique(y[-c(1:tab[1])])) != 1)
           stop("Incorrect coding for type='paired'. \n")
           diffxx <- xx2 - xx1
           low <- rowMin(diffxx)
           up  <- rowMax(diffxx)
           statisticlow <- ifelse(low > 0, low, 0)
           statisticup <- ifelse(up < 0, abs(up), 0)
           statistic <- pmax(statisticlow, statisticup)
          }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           low <- rowMin(x)
           up  <- rowMax(x)
           statisticlow <- ifelse(low > 0, low, 0)
           statisticup <- ifelse(up < 0, abs(up), 0)
           statistic <- pmax(statisticlow, statisticup)
          }

          pvals <- rep(NA, nrow(x))
          ranking <- order(abs(statistic), decreasing=TRUE)
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
           names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="Gapstatistic")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingGap", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL,...)
          RankingGap(x, y=as.numeric(y), type, gene.names, ...))

### signature: x=ExpressionSet, y=character.

setMethod("RankingGap", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL,...)
          RankingGap(exprs(x), y=pData(x)[,y], type, gene.names, ...))


