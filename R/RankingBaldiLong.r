### filename: RankingBaldiLong.r
### Title: Gene rankings via Baldi & Longs's Bayesian t statitic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 5.9.2007
### date(s) of updates: 6.9.2007
#
### Brief description:
#
#   Bayesian t-test following Baldi & Long (2001).
#   R code adaption from Anthony D. Long / Suman Sundaresh
#
#
### Further comments and notes:
#
#   s. also RankingFoxDimmic.r.
#   Produces results that differ slightly from the code above,
#   due to the sliding window problem. 
#
###**************************************************************************###

setGeneric("RankingBaldiLong", function(x, y, type=c("unpaired", "paired", "onesample"), m=100, conf=NULL,
            pvalues=TRUE, gene.names=NULL,...) standardGeneric("RankingBaldiLong"))


### signature: x=matrix, y=numeric.

setMethod("RankingBaldiLong", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
          m=100, conf=NULL, pvalues=TRUE, gene.names=NULL, ...){
          taby <- table(y)
          if(any(taby <= 2)) stop("More than two arrays for each level of y required. \n")
          ly <- length(y)
          if(ly != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          mode(x) <- "numeric"
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingBaldiLong)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample' \n")
          winsize <- ceiling(m/2)
          if(!is.null(conf)){
           conf <- as.integer(conf)
           if(conf < 0) stop("Argument 'conf' must be non-negative \n")
           }
          else{
                if(ly <= 4) conf <- 3*ly
                if(ly>4 & ly < 20) conf <- 2*ly
                else conf <- ly
              }
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
          m1 <- rowMeans(x)
          o1 <- order(m1)
          df1 <- ly-1
          ss1uo <- sqrt(rowSums((x-m1)^2)/df1)
          ss1 <- ss1uo[o1]
          p <- nrow(x)
          lower <- c(rep(1, winsize), 2:(p-winsize+1))
          upper <- c((winsize+1):p, rep(p, winsize))
          bayes <- sapply(1:p, function(i) {
                          low <- lower[i]
                          upp <- upper[i]
                          winlen <- length(low:upp)
                          ss11 <- sum(ss1[low:upp])/winlen
                          return(ss11)})
          part1 <- double(p)
          part1[o1] <- bayes
          bayesSD1 <- sqrt((conf * part1^2 + df1 * ss1uo^2)/(conf + df1-1))
          statistic <- m1/bayesSD1*sqrt(df1)
	        dft <- df1 + conf-1
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
           o1 <- order(m1)
           df1 <- (ly/2)-1
           ss1uo <- sqrt(rowSums((diffxx-m1)^2)/df1)
           ss1 <- ss1uo[o1]
           p <- nrow(x)
           lower <- c(rep(1, winsize), 2:(p-winsize+1))
           upper <- c((winsize+1):p, rep(p, winsize))
           bayes <- sapply(1:p, function(i) {
                          low <- lower[i]
                          upp <- upper[i]
                          winlen <- length(low:upp)
                          ss11 <- sum(ss1[low:upp])/winlen
                          return(ss11)})
           part1 <- double(p)
           part1[o1] <- bayes
           bayesSD1 <- sqrt((conf * part1^2 + df1 * ss1uo^2)/(conf + df1-1))
           statistic <- m1/bayesSD1*sqrt(df1)
	         dft <- df1+ conf-1
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
          o1 <- order(m1)
          o2 <- order(m2)
          df1 <- sum(ind)-1
          df2 <- sum(!ind)-1
          ss1uo <- sqrt(rowSums((xx1-m1)^2)/df1)
          ss2uo <- sqrt(rowSums((xx2-m2)^2)/df2)
          ss1 <- ss1uo[o1]
          ss2 <- ss2uo[o2]
          p <- nrow(x)
          lower <- c(rep(1, winsize), 2:(p-winsize+1))
          upper <- c((winsize+1):p, rep(p, winsize))
          bayes <- sapply(1:p, function(i) {
                          low <- lower[i]
                          upp <- upper[i]
                          winlen <- length(low:upp)
                          ss11 <- sum(ss1[low:upp])/winlen
                          ss22 <- sum(ss2[low:upp])/winlen
                          return(c(ss11, ss22))})
          part1 <- part2 <- double(p)
          part1[o1] <- bayes[1,]
          part2[o2] <- bayes[2,]
          bayesSD1 <- sqrt((conf * part1^2 + df1 * ss1uo^2)/(conf + df1-1))
		      bayesSD2 <- sqrt((conf * part2^2 + df2 * ss2uo^2)/(conf + df2 - 1))
          statistic <- (m1 - m2)/sqrt(((df1 * bayesSD1^2 + df2 * bayesSD2^2)/(df1 +
                df2)) * (((df1+1) + (df2+1))/((df1+1) * (df2+1))))
	        dft <- df1+df2 + 2*conf-2
          }

          if(pvalues) pvals <- 1 - pf(statistic^2, 1, dft)
          else pvals <- rep(NA, nrow(x))
          ranking <- order(abs(statistic), decreasing=TRUE)
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
          names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="BaldiLongT")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingBaldiLong", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"), m=50, conf=NULL,
            pvalues=TRUE, gene.names=NULL,...)
          RankingBaldiLong(x, y=as.numeric(y), type, m, conf, pvalues=TRUE, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingBaldiLong", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"), m=50, conf=NULL, pvalues=TRUE, gene.names=NULL,...)
          RankingBaldiLong(exprs(x), y=pData(x)[,y], type, m, conf, pvalues=TRUE, gene.names, ...)
          )
          
