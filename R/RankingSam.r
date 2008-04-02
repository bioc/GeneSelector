### filename: RankingSAM.r
### Title: Gene rankings via the Sam statistic of Tusher et al. (2001),
###        for a clear reference, prefer Schwender et al.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 6.9.2007
### date(s) of updates:
#
### Brief description:
#
#  Wrapper to 'samr'  from the package of the same name.
#
#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("RankingSam", function(x, y, type=c("unpaired", "paired", "onesample"),
            pvalues=TRUE, gene.names=NULL, ...) standardGeneric("RankingSam"))


### signature: x=matrix, y=numeric.

setMethod("RankingSam", signature(x="matrix", y="numeric"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
          pvalues=TRUE, gene.names=NULL, ...){
          require(samr, quietly=TRUE)
          mode(x) <- "numeric"
          if(length(y) != ncol(x)) 
          stop("Length of y is not equal to the number of columns of the expression matrix \n.")
          type <- match.arg(type)
          if( !is.element(type, eval(formals(RankingSam)$type)))
          stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
          yinput <- y
          taby <- table(y)
          if(is.element(type, c("unpaired", "paired")) & length(taby) != 2)
          stop("y has not exactly two levels ! \n")
          if( type == "paired"){
          if(taby[1] != taby[2] || length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
          stop("Incorrect coding for type='paired'. \n")
           y <- c(-(1:taby[1]),(1:taby[1]))
          }
          if(type == "onesample" & length(unique(y)) != 1){
          warning("Type has been chosen 'onesample', but y has more than one level. \n")
           y <- rep(1, length(y))
          }
          
          sam.type <- switch(type, unpaired ="Two class unpaired",
                                   paired = "Two class paired",
                                   onesample = "One class")
          if(type == "paired"){
          datobj <- list(x=x, y=y,
                         geneid=as.character(1:nrow(x)),
                         genenames=gene.names,
                         logged2=TRUE)
          }
          else{
          datobj <- list(x=x, y=as.numeric(as.factor(y)),
                         geneid=as.character(1:nrow(x)),
                         genenames=gene.names,
                         logged2=TRUE)
          }

          samtest <- samr(datobj, resp.type=sam.type, ...)
          statistic <- samtest$tt
          if(pvalues) pvals <- samr.pvalues.from.perms(statistic, samtest$ttstar)
          else pvals <- rep(NA, nrow(x))
          ranking <- order(abs(statistic), decreasing=TRUE)
          if(!is.null(gene.names))
            names(pvals) <- names(statistic) <- gene.names
          else{
          if(!is.null(rownames(x)))
            names(pvals) <- names(statistic) <- rownames(x)
          }
          new("GeneRanking", x=x, y=as.factor(yinput), statistic=statistic[ranking],
          ranking=ranking, pval=pvals[ranking], type=type, method="Sam")
})

### signature: x=matrix, y=factor.

setMethod("RankingSam", signature(x="matrix", y="factor"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                   pvalues=TRUE, gene.names=NULL, ...)
          RankingSam(x, y=as.numeric(y), type, pvalues, gene.names, ...)
          )

### signature: x=ExpressionSet, y=character.

setMethod("RankingSam", signature(x="ExpressionSet", y="character"),
          function(x, y, type=c("unpaired", "paired", "onesample"),
                   pvalues=TRUE, gene.names=NULL, ...)
          RankingSam(exprs(x), y=pData(x)[,y], type, pvalues, gene.names, ...)
          )

