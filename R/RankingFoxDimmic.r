### filename: RankingFoxDimmic.r
### Title: Gene rankings via Fox & Dimmic's Bayesian t statitic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 4.9.2007
#
### Brief description:
#
#   Bayesian t-test following Fox & Dimmic (2006).
#   R code adaption from Richard J. Fox Code
#
#
### Further comments and notes:
#
#   Produces exactly the same results as the original code.
#   s. also RankingBaldiLong.r.
#
###**************************************************************************###

setGeneric("RankingFoxDimmic", function(x, y, type="unpaired", m=4,
            pvalues=TRUE, gene.names=NULL,...) standardGeneric("RankingFoxDimmic"))


### signature: x=matrix, y=numeric.

setMethod("RankingFoxDimmic", signature(x="matrix", y="numeric"),
          function(x, y, type="unpaired", m=4, pvalues=TRUE, gene.names=NULL, ...){
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingFoxDimmic)$type)))
          stop("Currently, only type='unpaired' is implemented. \n")
          win.upper <- ceiling(m/2)
          win.lower <- ceiling(-m/2)
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
          o1 <- order(m1)
          o2 <- order(m2)
          ss1 <- rowSums((xx1-m1)^2)[o1]
          ss2 <- rowSums((xx2-m2)^2)[o2]
          df1 <- sum(ind)-1
          df2 <- sum(!ind)-1
          p <- nrow(x)
          lower <- (1:p)+win.lower
          upper <- (1:p)+win.upper
          indlow <- (lower<1)
          lower[indlow] <- 1
          upper[indlow] <- m+1
          indup <- (upper>p)
          lower[indup] <- p-m
          upper[indup] <- p
          bayes <- sapply(1:p, function(i) {
                          low <- lower[i]
                          upp <- upper[i]
                          ss11 <- sum(ss1[low:upp])
                          ss22 <- sum(ss2[low:upp])
                          return(c(ss11, ss22))})
          part1 <- part2 <- double(p)
          part1[o1] <- bayes[1,]
          part2[o2] <- bayes[2,]
          ss.bayes <- part1 + part2
          df.bayes <- (m+1)*(df1 + df2)
          sd.bayes <- sqrt(ss.bayes/df.bayes)
          delta <- m1 - m2
          statistic <- delta/(sd.bayes*sqrt(1/(df1+1) + 1/(df2+1)))
          }

          if(pvalues) pvals <- 1 - pf(statistic^2, 1, df.bayes)
          else pvals <- rep(NA, nrow(x))
          ranking <- rank(-abs(statistic), ties = "first")
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
          names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y = as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="FoxDimmicT")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingFoxDimmic", signature(x="matrix", y="factor"),
          function(x, y, type="unpaired", m=4,
            pvalues=TRUE, gene.names=NULL,...)
          RankingFoxDimmic(x, y=as.numeric(y), type, m, pvalues=TRUE, gene.names, ...)

          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingFoxDimmic", signature(x="ExpressionSet", y="character"),
          function(x, y, type="unpaired", m=4,
            pvalues=TRUE, gene.names=NULL,...)
          RankingFoxDimmic(exprs(x), y=pData(x)[,y], type, m, pvalues=TRUE, gene.names, ...)

          )


