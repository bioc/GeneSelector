### filename: RankingTstat.r
### Title: Gene rankings via the t-statistic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.8.2007
### date(s) of updates: 24.8.2007, 28.8.2007, 10.9.2007
#
### Brief description:
#
#   One of a lot similarly structured ranking methods.
#
### Further comments and notes:
#
#   24.8.: correction. P-values are now ordered !
#   28.8.: computation of pvalues now optional.
#   10.9.: Complete new version in order to remove gene filter dependency. 
#
###**************************************************************************###


### generic function

setGeneric("RankingTstat", function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL,...)
            standardGeneric("RankingTstat"))
            
### signature: x=matrix, y=numeric.

setMethod("RankingTstat", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                  pvalues=TRUE, gene.names=NULL, ...){
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingTstat)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
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
           ss2 <- rowSums((xx2-m2)^2)
           s <- sqrt((ss1+ss2)/(df1+df2))
           se <- s*sqrt(1/(df1+1)+ 1/(df2+1))
           ttest <- (m1-m2)/se
           if(pvalues) pvals <- 1-pf(ttest^2, 1, df1+df2)
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
           m <- rowMeans(diffxx)
           ss <- rowSums((diffxx-m)^2)
           df <- length(y)/2-1
           ttest <- m/sqrt(ss/df)*sqrt(df+1)
           if(pvalues) pvals <- 1 - pf(ttest^2, 1, df)
           else pvals <- rep(NA, nrow(x))
          }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
          m <- rowMeans(x)
          ss <- rowSums((x-m)^2)
          df <- length(y)-1
          ttest <- m/sqrt(ss/df)*sqrt(df+1)
          if(pvalues) pvals <- 1 - pf(ttest^2, 1, df)
          else pvals <- rep(NA, nrow(x))
          }
          statistic <- ttest
          ranking <- rank(-abs(statistic))
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y = as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type,
          method="ordinaryT")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingTstat", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL,...)
          RankingTstat(x, y=as.numeric(y), type, pvalues=TRUE, gene.names, ...))

### signature: x=ExpressionSet, y=character.

setMethod("RankingTstat", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL,...)
          RankingTstat(exprs(x), y=pData(x)[,y], type,
          pvalues=TRUE, gene.names, ...))

