### filename: RankingEbam.r
### Title: Gene rankings via Efron's empirical bayes mixture model approach
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.8.2007
#
### Brief description:
#
#   A wrapper to the siggenes package.
#
#
#
### Further comments and notes:
#
#   s. also RankingWilcEbam.r for a rank-based version of this test.
#
###**************************************************************************###

### library(Biobase) ->>> Definition of class 'ExpressionSet'

### generic function

setGeneric("RankingEbam", function(x, y, type=c("unpaired", "paired", "onesample"),
            gene.names=NULL,...) standardGeneric("RankingEbam"))


### signature: x=matrix, y=numeric.

setMethod("RankingEbam", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                  gene.names=NULL, p0.estimation = c("splines", "interval", "adhoc"), ...){
          require(siggenes, quietly=TRUE)
          ll <- eval(substitute(list(...)))
          if(!hasArg(quan.a0)) ll$quan.a0 <- 0.9
          if(!hasArg(p0)) p0 <- NA
          if(!hasArg(lambda)) lambda <- NULL
          if(!hasArg(use.weights)) use.weights <- FALSE
          if(!hasArg(ncs.value)) ncs.value <- "max"
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          ll$data <- x
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingEbam)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
          y <- as.factor(y)
          if(type == "unpaired"){
          if(nlevels(y) != 2)
          stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
           ll$cl <- as.numeric(y)
           out <- do.call(z.ebam, ll)
          }
          if(type == "paired"){
           tab <- table(y)
           if(length(tab) != 2)
           stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
           if(tab[1] != tab[2] || length(unique(y[1:tab[1]])) != 1 | length(unique(y[-c(1:tab[1])])) != 1)
           stop("Incorrect coding for type='paired'. \n")
           ll$cl <- c(-(1:tab[1]),(1:tab[1]))
           out <- do.call(z.ebam, ll)
          }
          
          if(type == "onesample"){
          if(length(unique(y)) != 1)
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
          ll$cl <- rep(1, length(y))
          out <- do.call(z.ebam, ll)
          }
          
          check.out <- siggenes:::checkFUNout(out)
          if (is.na(p0))
          p0 <- siggenes:::pi0.est3(out, p0.estimation, exact = check.out$exact,
                          lambda = lambda, ncs.value = ncs.value, use.weights = use.weights)
          if (!is.null(gene.names)) {
          names(out$z) <- names(out$ratio) <- substring(gene.names, 1, 50)
          if (length(check.out$vec.pos) != 0)
              names(check.out$vec.pos) <- names(check.out$vec.neg) <- names(out$z)
          }
        posterior <- 1 - p0 * out$ratio
        posterior[posterior < 0] <- 0
        statistic <- posterior
        ranking <- order(posterior, decreasing=TRUE)
        pvals <- rep(NA, length(posterior))
          if(!is.null(gene.names))
          names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
          names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(y), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="Ebam")
}
)

### signature: x=matrix, y=factor.

setMethod("RankingEbam", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
             gene.names=NULL, p0.estimation = c("splines", "interval", "adhoc"), ...)
          RankingEbam(x, y=as.numeric(y), type, gene.names, p0.estimation, ...)
        )

### signature: x=ExpressionSet, y=character.

setMethod("RankingEbam", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
             gene.names=NULL, p0.estimation = c("splines", "interval", "adhoc"), ...)
          RankingEbam(exprs(x), y=pData(x)[,y], type, gene.names, p0.estimation,...))


