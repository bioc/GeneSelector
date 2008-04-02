### filename: RankingPermuationT.r
### Title: Gene rankings via permuation t stastic.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### date(s) of updates:
#
### Brief description:
#
#  Wrapper to "mt.sample.teststat"  from the package "multtest"
#  that uses a fast C routine. Exclusively possible for
#  type = "unpaired" (two sample case).
#
#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("RankingPermutation", function(x, y, type="unpaired", B=100,
            gene.names=NULL, ...) standardGeneric("RankingPermutation"))


### signature: x=matrix, y=numeric.

setMethod("RankingPermutation", signature(x="matrix", y="numeric"),
          function(x, y, type="unpaired", B=100, gene.names=NULL, ...){
          require(multtest, quietly=TRUE)
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingPermutation)$type)))
          stop("Permutation test is only possible for the two class
                unpaired setting \n")
          ll <- eval(substitute(list(...)))
          taby <- table(y)
          if(length(taby) != 2)
          stop("y has not exactly two levels ! \n")
          if(!hasArg(test)) ll$test <- "t.equalvar"
          if(B>1000){
            B <- 1000
            warning("number of permutations > 1000; reset to 1000 \n")
            }
          ll$B <- B
          ll$classlabel <- ifelse(y==names(taby[2]),0,1)
          permute <- apply(x, 1, function(z){
                        ll$V <- z
                        res <- do.call(mt.sample.teststat, ll)
                        statistic <- res[1]
                        perm <- res[-1]
                        pvalue <- mean(abs(statistic)<abs(perm))
                        return(c(statistic, pvalue))
          })
          statistic <- permute[1,]
          pvals <- permute[2,]

          ranking <- order(pvals, 1/abs(statistic))
          if(!is.null(gene.names))
            names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="Permutation")
})

### signature: x=matrix, y=factor.

setMethod("RankingPermutation", signature(x="matrix", y="factor"),
          function(x, y, type="unpaired", B=100, gene.names=NULL, ...)
          RankingPermutation(x, y=as.numeric(y), type, B, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingPermutation", signature(x="ExpressionSet", y="character"),
          function(x, y, type="unpaired", B=100, gene.names=NULL, ...)
          RankingPermutation(exprs(x), y=pData(x)[,y], type, B, gene.names, ...)
          )

