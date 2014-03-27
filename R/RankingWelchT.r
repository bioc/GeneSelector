### filename: RankingWelchT.r
### Title: Gene rankings via the Welch t-statistic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation:  10.9.2007
### date(s) of updates:
#
### Brief description:
#
#   Welch's t-Test for two-sample unpaired data.
#
### Further comments and notes:
#
#   Originally included in RankingTstat, for compatibility issues removed
#   and defined as independent method.
#
###**************************************************************************###


### generic function

setGeneric("RankingWelchT", function(x, y, type="unpaired",
            pvalues=TRUE, gene.names=NULL,...) standardGeneric("RankingWelchT"))

### signature: x=matrix, y=numeric.

setMethod("RankingWelchT", signature(x="matrix", y="numeric"),
          function(x, y, type="unpaired",
                  pvalues=TRUE, gene.names=NULL, ...){
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingWelchT)$type)))
          stop("'type' must be 'unpaired'. \n")
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
           se1 <- ss1/(df1*(df1+1))
           se2 <- ss2/(df2*(df2+1))
           se12 <-  sqrt(se1 + se2)
           ttest <- (m1-m2)/se12
           term <- se1/(se1+se2)
           nu <- floor(1/((term^2/df1 + (1-term)^2/df2)))
           if(pvalues) pvals <- 1 - pf(ttest^2, 1, nu)
           else pvals <- rep(NA, nrow(x))
          statistic <- ttest
          ranking <- rank(-abs(statistic), ties = "first")
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y = as.factor(y), statistic=statistic,
          ranking=ranking, pval=pvals, type=type, method="WelchT")})

### signature: x=matrix, y=factor.

setMethod("RankingWelchT", signature(x="matrix", y="factor"),
          function(x, y, type="unpaired",
            pvalues=TRUE, gene.names=NULL,...)
          RankingWelchT(x, y=as.numeric(y), type, pvalues=TRUE, gene.names, ...))

### signature: x=ExpressionSet, y=character.

setMethod("RankingWelchT", signature(x="ExpressionSet", y="character"),
          function(x, y, type="unpaired",
            pvalues=TRUE, gene.names=NULL,...)
          RankingWelchT(exprs(x), y=pData(x)[,y], type,
          pvalues=TRUE, gene.names, ...))
          
