### filename: RankingWilcoxon.r
### Title: Gene rankings via the Wilcoxon statistic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.8.2007
#
### Brief description:
#
#   Ranking based on the function ranks() from R
#   instead of using the slow Wilcoxon-Test
#   of p-values possible. ~ 1.10-15 seconds one execution
#
### Further comments and notes:
#
#
#
###**************************************************************************###

### library(Biobase) ->>> Definition of class 'ExpressionSet'

### generic function

setGeneric("RankingWilcoxon", function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=FALSE, gene.names=NULL,...) standardGeneric("RankingWilcoxon"))


### signature: x=matrix, y=numeric.

setMethod("RankingWilcoxon", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=FALSE, gene.names=NULL, ...){
          mode(x) <- "numeric"
          if(length(y) != ncol(x))
          stop("Length of y is not equal to the number of columns of the expression matrix \n.") 
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingWilcoxon)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
          y <- as.factor(y)
          if(type == "unpaired"){
          if(nlevels(y) != 2)
          stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
          taby <- table(y)
          levy <- names(taby)[which.max(taby)]
          ind <- (y==levy)
          Rx <- apply(x, 1, rank)
          r1 <- colSums(Rx[ind, ,drop=FALSE])-sum(1:max(taby))
          e1 <- taby[1]*taby[2]/ 2
          #browser()                       
          if(pvalues){
              maxr <- sum((min(taby)+1):length(y))-sum(1:sum(ind)) 
              pvals <- 2*pwilcox(ifelse(r1<e1, maxr-r1, r1), 
                                taby[1], taby[2], lower.tail=FALSE)
            }
                                 
          else pvals <- rep(NA, nrow(x))
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
           #indx <- apply(diffxx, 1, ">", 0)
           r1 <- apply(diffxx, 1, function(z){zz <- rank(abs(z)); sum(zz[z>0])})
           ly <- length(y)
           #r1 <- colSums(Rx[diffx>0, ,drop=FALSE])
           e1 <- (ly/2)*(ly/2+1)/4
           if(pvalues){
             maxr <- sum(1:(ly/2))
             pvals <- 2*psignrank(ifelse(r1<e1, maxr-r1, r1), n=ly/2, lower.tail=FALSE)
             }
           else pvals <- rep(NA, nrow(x))
          }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
          #indx <- apply(x, 1, ">", 0)
          r1 <- apply(x, 1, function(z){ zz <- rank(abs(z)); sum(zz[z>0])})
          ly <- length(y)
          #r1 <- colSums(Rx[indx, ,drop=FALSE])
          e1 <- (ly)*(ly+1)/4
          maxr <- sum(1:ly)                        
          if(pvalues){
                pvals <- 2*psignrank(ifelse(r1<e1, maxr-r1, r1), n=ly, lower.tail=FALSE)
            }
           else pvals <- rep(NA, nrow(x))
          }

          statistic <- r1
          ranking <- rank(-abs(r1 - e1), ties.method = "first")
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
          names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y = as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="Wilcoxon")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingWilcoxon", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=FALSE, gene.names=NULL,...)
          RankingWilcoxon(x, y=as.numeric(y), type, pvalues, gene.names,...)

          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingWilcoxon", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=FALSE, gene.names=NULL,...)
          RankingWilcoxon(exprs(x), y=pData(x)[,y], type, pvalues, gene.names,...)

          )


