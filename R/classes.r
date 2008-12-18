###**************************************************************************###
### filename: classes.r
### Title: Raw versions of classes.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.8.2007
### date(s) of updates:
#   17.8.2007; 23.8.2007 ; 27.8.2007 ; 28.8.2007; 3.9.2007; 10.9.2007; 11.9.2007;
#   25.9.2007; 26.9.2007
### Brief description:
#
#   Stores the class definitions and convenience methods for them,
#   is frequently updated.
#
### Further comments and notes:
#
#
#
###**************************************************************************###

#+++++++++++ Class: GeneRanking +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))

### show-method: can still be improved w.r.t to method output
### -> will be kept for simplicity -- would require 'switch'
                        
setMethod("show", signature="GeneRanking", function(object){
          cat("Ranking by ", object@method, "\n", sep="")
          cat("number of genes: ", length(object@statistic), "\n", sep="")
          })

### summary-method

setGeneric("summary", function(object, ...) standardGeneric("summary"))
          
setMethod("summary", signature="GeneRanking", function(object){
          if(is.na(object@pval[1])){ 
           summ <- as.matrix(summary(object@statistic))
           colnames(summ) <- "statistic"
           print(summ)
           }
          else cbind(statistic=summary(object@statistic), p_values=summary(object@pval))}
          )
          
### toplist-method

setGeneric("toplist", function(object, top=10, show=TRUE,...) standardGeneric("toplist"))

setMethod("toplist", signature(object="GeneRanking"), function(object, show=TRUE, top=10){
          ind <- object@ranking[1:top]
          ret <- data.frame(index=ind, statistic=object@statistic[1:top], pvals=object@pval[1:top])
          rownames(ret) <- NULL
          if(show) print(ret)
          invisible(ret)
          })
          
#+++++++++++ Class: FoldMatrix ++++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="FoldMatrix",
        representation(foldmatrix="matrix", k="integer", replicates="integer",
                       type="character", minclassize="integer", balanced="logical"))
                       
### show-method:

setMethod("show", signature="FoldMatrix", function(object){
          if(object@type == "paired")
          cat("number of removed samples per replicate: ", 2*object@k, "\n", sep="")
          else
          cat("number of removed samples per replicate: ", object@k, "\n", sep="")
          cat("number of replicates: ", object@replicates, "\n", sep="")
          if(object@type=="unpaired"){
           if(object@balanced) cat("constraints:  balanced subsampling. \n")
           else cat("constraints: minimum classize for each class: ", object@minclassize, "\n", sep="")
           }
          })
          
setMethod("summary", signature="FoldMatrix", function(object, repl=1){
          FF <-  object@foldmatrix
          if(any(repl > ncol(FF)))
          stop("Invalid value for 'repl' \n")
          cat("\t \t observations left out \n")
          for(i in seq(along=repl)){
          cat("replication ", repl[i], " : ")
          whichout <- which(FF[,repl[i]] == FALSE)
          for(j in 1:length(whichout)){
          cat(whichout[j], ", ")
          }
          cat("\n")}})
          
#+++++++++++ Class: BootMatrix ++++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="BootMatrix",
        representation(bootmatrix="matrix", replicates="integer",
                       type="character", maxties="vector", minclassize="integer", 
                        balancedclass="logical", balancedsample="logical"))
                       
### show-method:

setMethod("show", signature="BootMatrix", function(object){
          cat("number of bootstrap replicates: ", object@replicates, "\n", sep="")
          if(object@type=="unpaired"){
           if(object@balancedclass){
            if(object@balancedsample){
             if(is.na(object@maxties))
             cat("constraints:  balanced class and balanced sample resampling. \n")
             else{
              cat("constraints: balanced class and balanced sample resampling, 
                  maximum number of ties per observation: ",
                  object@maxties, "\n", sep="")}
               }
              else{
             if(is.na(object@maxties))
             cat("constraints:  balanced class resampling. \n")
             else{
              cat("constraints: balanced class resampling, 
                  maximum number of ties per observation: ",
                object@maxties, "\n", sep="")}
               } 
             }
              else{
                if(object@balancedsample){ 
                if(is.na(object@maxties))
                cat("constraints: balanced sample resampling,
                     minimum classize for each class: ", object@minclassize, "\n", sep="")
                else{
                cat("constraints: balanced sample resampling, minimum classize for each class: ", 
                    object@minclassize, "\n \t     maximum number of ties per observation: ",
                 object@maxties, "\n", sep="")}
              }
              else{
               if(is.na(object@maxties))
                cat("constraints: minimum classize for each class: ", 
                     object@minclassize, "\n", sep="")
                else{
                cat("constraints: minimum classize for each class: ", 
                     object@minclassize, "\n \t     maximum number of ties per observation: ",
                      object@maxties, "\n", sep="")} 
              }
             }
            }
            else{
             if(object@balancedsample){
             if(is.na(object@maxties))
              cat("constraints: balanced sample resampling \n", sep="")
             else
              cat("constraints: balanced sample resampling, maximum number of ties per observation: ", object@maxties, "\n", sep="")
             }
             else{
             if(!is.na(object@maxties))
              cat("maximum number of ties per observation: ", object@maxties, "\n", sep="")
            }
           }
        } 
       )
       
setMethod("summary", signature="BootMatrix", function(object, repl=1){
          BB <-  object@bootmatrix
          if(any(repl > ncol(BB)))
          stop("Invalid value for 'repl' \n")
          tabmatrix <- t(apply(BB, 2, function(z) tabulate(z, nrow(BB))))
          rownames(tabmatrix) <- paste("repl.", 1:ncol(BB), sep="")
          colnames(tabmatrix) <- paste("#obs.", 1:nrow(BB), sep="")
          print(tabmatrix[repl,])
          })
              
#+++++++++++ Class: RepeatRanking +++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="RepeatRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))
                       
setMethod("show", signature="RepeatRanking", function(object){
          cat(ncol(object@rankings), " rankings with scheme '", object@scheme, "' \n", sep="")
          cat("Method used: ", object@original@method, "\n", sep="")
          })
          
setMethod("toplist", signature(object="RepeatRanking"), function(object, top=10, show=TRUE){
          inds <- apply(object@rankings, 2, function(z) z[1:top])
          uniqinds <- unique(as.vector(inds))
          emptytab <- numeric(length(uniqinds))
          names(emptytab) <- as.character(uniqinds)
          factinds <- apply(inds, 1, factor, levels=uniqinds)
          genetable <- apply(factinds, 2, function(z){ tab <- table(z) 
                                                       emptytab[names(tab)] <- tab[names(tab)]
                                                       return(emptytab)})
          rown <- rownames(genetable)
          ord <- order(genetable[,1], decreasing=TRUE)
          genetableord <- genetable[ord,]
          printgenetable <- data.frame(genetableord)
          rownames(genetableord) <- rown
          colnames(printgenetable) <- paste("Rank ", 1:top, sep="")     
          cat("original dataset: \n \n")
          orig <- object@original@ranking[1:top]
          print(data.frame(index=orig, statistic=object@original@statistic[1:top], pvals=object@original@pval[1:top]))
          cat("\n \n")
          toplists <- vector(mode="list", length=ncol(object@rankings))
           
          for(i in 1:ncol(object@rankings)){
          toplists[[i]] <- data.frame(index=inds[,i], statistic=object@statistics[1:top,i], 
                                      pvals=object@pvals[1:top,i])
          if(show){
              cat("Replication", i, ": \n \n")
              print(toplists[[i]])
              cat("\n \n")
            }
          }
          if(show){ 
          cat("In the following table, rownames correspond to gene indices. \n") 
          cat("The columns contain the absolute frequencies for the corresponding ranks \n")
          cat("over all replications. \n") 
          cat("Genes are ordered according to the first column, then to the second, and so on. \n \n")
          print(printgenetable)
          cat("\n \n")
          }
          return(invisible(toplists))})
          
setGeneric("plot")

setMethod("plot", signature("RepeatRanking", "missing"),
          function(x, frac=1/100, ...){
          if(frac < 0 | frac > 1) stop("frac must be between 0 and 1. \n")
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          if(!hasArg(xlab)) ll$xlab <- "Ranks in the orginal dataset"
          if(!hasArg(ylab)) ll$ylab <- "Ranks in perturbed datasets"
          if(!hasArg(cex))  ll$cex <- 0.5
          if(!hasArg(main)) ll$main <- ""
          ngenes <- ceiling(nrow(x@rankings)*frac)
          ll$x <- ll$y <- 1:ngenes 
          ll$type <- "n"
          ll$ylim <- c(0, 2*ngenes)
          do.call(plot, args=ll)
          abline(0,1, lwd=2.5, col="blue")
          ll$type <- "p"
          for(i in 1:ncol(x@rankings)){
          ll$y <- match(x@original@ranking[1:ngenes], 
                       x@rankings[,i])
          do.call(points, args=ll)               
          }
          }
          )
          
### join-Method still needs testing ! ->> done

setGeneric("join", function(RR1, RR2) standardGeneric("join"))
                  
setMethod("join", signature("RepeatRanking", "RepeatRanking"),
           function(RR1, RR2){
          r1 <- RR1@rankings
          r2 <- RR2@rankings
          if(!identical(RR1@original, RR2@original))
          stop("Input arguments do not match, different original ranking \n")
          if(nrow(r1) != nrow(r2)) 
          stop("Input arguments do not match, different number of rows in ranking matrix \n")
          if(RR1@original@method != RR2@original@method)
          stop("Input arguments do not match, different ranking methods used")
           r12 <- cbind(r1, r2)
           t12 <- cbind(RR1@statistics, RR2@statistics)
           if(is.na(RR1@pvals)[1]) p12 <- NA 
           else p12 <- cbind(RR1@pvals, RR2@pvals)
           lr1 <- ncol(r1)
           lr2 <- ncol(r2)
           colnames(r12) <- colnames(p12) <- colnames(t12) <- rep("", (lr1+lr2))
           colnames(r12)[1:lr1] <- colnames(p12)[1:lr1] <- colnames(t12)[1:lr1] <- rep(RR1@scheme, lr1)
           colnames(r12)[-(1:lr1)] <- colnames(p12)[-(1:lr1)] <- colnames(t12)[-(1:lr1)] <- rep(RR2@scheme, lr2) 
           scheme <- "combined"
           new("RepeatRanking", original=RR1@original, rankings=r12, pvals=p12, statistics = t12, scheme=scheme)
          } 
         )
         
         
### needs still documentation !!!

setGeneric("variance", function(RR, estimator = c("var", "mad", "iqr"), center = c("perturbed", "original"))
            standardGeneric("variance"))

setMethod("variance", signature("RepeatRanking"),
           function(RR, estimator = c("var", "mad", "iqr"), center = c("perturbed", "original")){
           ranking <- RR@rankings
           ranking <- apply(ranking, 2, function(z) match(1:nrow(ranking), z))
           if(ncol(ranking) < 2) stop("Too few replications in the argument 'RR'; variance cannot be estimated \n")
           ord <- r0 <- RR@original@ranking
           r0 <- match(1:nrow(ranking), r0)
           estimator <- match.arg(estimator)
           center <- match.arg(center)
           if(!is.element(estimator, c("var", "mad", "iqr")))
           stop("Invalid 'estimator' specified \n")
           if(!is.element(center, c("perturbed", "original")))
           stop("Invalid 'center' specified \n")
           if(estimator == "iqr"){
            iqrs <- apply(ranking, 1, IQR)[ord]
            names(iqrs) <- paste("top gene", 1:nrow(ranking)) 
            return(iqrs)
           }
            if(estimator == "var"){
             if(center == "perturbed")
              vars <- apply(ranking, 1, var)[ord]
             else{
              vars <- 1/(ncol(ranking)-1)*rowSums((ranking - r0)^2)[ord]
              }
              names(vars) <- paste("top gene", 1:nrow(ranking)) 
              return(vars)
            }
            if(estimator == "mad"){
             if(center == "perturbed")
              mads <- ((apply(ranking, 1, mad, constant = 1))^2)[ord]
             else{
              mads <- ((apply(abs(ranking - r0), 1, median))^2)[ord]
             }
             names(mads) <- paste("top gene", 1:nrow(ranking)) 
             return(mads)
          }
          })

          
#+++++++++++ Class: StabilityLm +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="StabilityLm",
        representation(coefficients="numeric", R2vec="numeric",
                       multivariateR2 = "numeric", residuals = "numeric",
                       residuals.unscaled = "numeric", weightscheme = "list"))
                       
setMethod("show", signature="StabilityLm", function(object){
          wl <- object@weightscheme
          if(wl$scheme == "pval") scheme <- "based on p-values"
          else scheme <- "based on ranks"
          decay <- switch(wl$decay, linear="linear weight decay",
                                    quadratic="quadratic weight decay",
                          exponential=paste("exponential weight decay, alpha=", wl$alpha))
          
          cat("Stability measure: weighted linear regression, \n",
               "weighting: ", scheme, ", ", decay, "\n",
               "multivariate R2 is: ", object@multivariateR2, "\n") 
          })
          
setMethod("plot", signature("StabilityLm", "missing"), 
           function(x, frac=1/50, scaled = TRUE, standardize=TRUE,...){
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          if(!hasArg(xlab)) ll$xlab <- "Ranks of the orginal dataset"
          if(!hasArg(ylab)) ll$ylab <- "Multivariate Regression residual"
          if(!hasArg(cex))  ll$cex <- 0.8
          if(!hasArg(pch)) ll$pch <- "o"
          if(!hasArg(main)) ll$main <- ""
          ngenes <- ceiling(length(x@residuals)*frac)
          ll$x <- 1:ngenes
          if(scaled) ll$y <- x@residuals[1:ngenes]
          else ll$y <- x@residuals.unscaled[1:ngenes]
          if(standardize) ll$y <- scale(ll$y)
          do.call(plot, args=ll)
          })
          
#+++++++++++ Class: StabilityOverlap +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="StabilityOverlap",
        representation(overlap = "matrix", scores ="matrix", weightscheme = "list"))
        
setMethod("show", signature="StabilityOverlap", function(object){
          wl <- object@weightscheme
          if(wl$scheme == "pval") scheme <- "based on p-values"
          else scheme <- "based on ranks"
          decay <- switch(wl$decay, linear="linear weight decay",
                                    quadratic="quadratic weight decay",
                          exponential=paste("exponential weight decay, alpha=", wl$alpha))
          p <- nrow(object@scores)
          cat("Stability measure: Overlap Score, \n",
               "weighting: ", scheme, ", ", decay, "\n",
               "mean Overlap Score is: ", mean(object@scores[p,]), "\n") 
          })
          
### for generic method definition, s. above
          
setMethod("summary", signature="StabilityOverlap", 
          function(object, which=c("scores", "overlap"), position=100,...){
          which <- match.arg(which)
          if(!is.element(which, c("scores", "overlap")))
          stop("'which' must be either 'scores' or 'overlap' \n") 
          position <- min(position, nrow(object@scores))
          if(which == "scores"){
            cat("summary of scores, iterationwise, up to rank", position,  ": \n \n")
            summary(object@scores[position,])
            }
          else {
           cat("summary of overlaps, iterationwise, up to rank", position, ": \n \n")
            summary(object@overlap[position,])
           }
           })
          
setMethod("plot", signature("StabilityOverlap", "missing"), 
           function(x, frac=1/50, ...){
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          ll$xlab <- "Rank"
          ll$ylab <- "Overlap"
          ll$lwd <- 1.2
          ll$main <- "Mean Overlap"
          ll$type <- "h"
          ll$col <- "grey"
          ngenes <- ceiling(nrow(x@scores)*frac)
          ll$x <- 1:ngenes
          ll$y <- apply(x@overlap[1:ngenes, , drop=FALSE], 1, mean)
          ll$ylim <- c(0, ngenes)
          layout(mat=as.matrix(c(1,2)), heights=c(1,1))
          do.call(plot, args=ll)
          abline(0,1, lwd=2.5)
          ll$col <- "red"
          ll$lty <- "dashed"
          ll$lwd <- 1.2
          ll$type <- "l"
          ll$y <- apply(x@scores[1:ngenes, , drop=FALSE], 1, mean)
          ll$ylim <- c(0,1)
          ll$ylab <- "score"
          ll$main <- "Mean Score"
          do.call(plot, args=ll)
          layout(mat=as.matrix(1))
           }
          )
          
         
#+++++++++++ Class: AggregatedRanking +++++++++++++++++++++++++++++++++++++++++#

setClass(Class="AggregatedRanking",
        representation(posterior = "vector", summary = "numeric", pval = "vector", type= "character",
                       fun="character", method = "character"))
        
           
setMethod("show", signature="AggregatedRanking", function(object){
            cat(object@type, " aggregation for ", length(object@summary), "genes \n")
            cat("aggregation statistic: ",  object@fun, "\n")
            })
            
setMethod("plot", signature("AggregatedRanking", "missing"), 
           function(x, index, ...){
          if(x@type == "simple") 
           stop("Plot only possible for bayesian aggregation \n")
          post <- x@posterior[[index]] 
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          ll$x <- as.numeric(names(post))
          ll$y <- post
          if(!hasArg(lwd)) ll$lwd <- 1.5
          if(!hasArg(main)) ll$main <- paste("posterior distribution for gene index ", 
                                            substitute(index))
          ll$type <- "h"
          if(!hasArg(ylab)) ll$ylab <- "posterior probability"
          if(!hasArg(xlab)) ll$xlab <- "rank"
          do.call(plot, args=ll)
           }
          )
          
#+++++++++++ Class: CombinedRanking +++++++++++++++++++++++++++++++++++++++++#

setClass(Class="CombinedRanking",
        representation(ranking = "numeric", rankmatrix = "matrix", inout = "matrix", 
                       selected = "numeric", adjpval="vector",  maxrank = "numeric", 
                       statistics = "character", absdist = "numeric", reldist = "numeric"))
        
           
setMethod("show", signature="CombinedRanking", function(object){
            cat("GeneSelector run with gene rankings from the following statistics: \n")
            for(i in 1:length(object@statistics)) cat(object@statistics[i], "\n")
            cat("Number of genes below threshold rank ", 
                 object@maxrank, " in all statistics:", 
                 length(object@selected), "\n")})
                 
setMethod("toplist", signature=(object="CombinedRanking"), function(object, show = TRUE, top=10){
          ind <- object@ranking[1:top]
          ret <- data.frame(index=ind, pvals=object@adjpval[1:top])
          rownames(ret) <- NULL
          if(show) print(ret)
          invisible(ret)
          })
          
setGeneric("SelectedGenes", function(object) standardGeneric("SelectedGenes"))

setMethod("SelectedGenes", signature(object="CombinedRanking"), function(object){
          toplist(object, length(object@selected)) })
          
          
setMethod("plot", signature("CombinedRanking", "missing"), 
           function(x, top=10, ...){
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          ll$height <- x@reldist[1:top]
          deltay <- max(ll$height)/top
          ll$horiz <- FALSE
          if(!hasArg(main)) ll$main <- paste("GeneSelector relative distance plot")
          if(!hasArg(ylab)) ll$ylab <- "relative distance"
          if(!hasArg(xlab)) ll$xlab <- "gene index"
          if(!hasArg(cex.lab)) ll$cex.lab <- 1.5
          if(!hasArg(ylim)) ll$ylim <- c(0, max(ll$height)+deltay*2)
          bb <- do.call(barplot, args=ll)
          chars <- as.character(x@ranking[1:top])
          for(i in 1:top) 
          characterplot(chars[i], bb[i], ll$height[i], deltax=3/top, deltay=deltay, cex=15/top)} 
          )


setGeneric("GeneInfoScreen", function(object, which) standardGeneric("GeneInfoScreen"))          

setMethod("GeneInfoScreen", signature(object="CombinedRanking"), 
          function(object, which){
          if(length(which) != 1) stop("Length of 'which' must be one. \n")
          stats <- rev(object@statistics)
          nostats <- length(stats)
          inout <- rev(object@inout[which,])
          RR <- rev(object@rankmatrix[which,])
          p <- nrow(object@rankmatrix)
          p0 <- c(0.01, 0.05, 0.1, 0.25, 1)
          defbreaks <- ceiling(p0*p)
          plot(0:1, 0:1, axes=FALSE, type="n", 
              main=paste("GeneInfoScreen for gene ", substitute(which)),
              xlab="", ylab="")
          partition <- seq(from=0, to=0.9, length=nostats+1)
          cols <- terrain.colors(5)
          for(i in 1:nostats){
          rect(0, partition[i], 0.25, partition[i+1])
          points(0.125, (partition[i] + partition[i+1])/2, cex=1.7, pch=inout)
          text(0.15, 1, "selected ?", cex=1)
          text(0.5, 1,  "statistic", cex=1)
          text(0.8, 1, "rank", cex=1)
          text(0.5, (partition[i] + partition[i+1])/2, stats[i], cex=1.5)
          text(0.8, (partition[i] + partition[i+1])/2, round(RR[i]), cex=1.5,
               col = cols[which(defbreaks >= RR[i])[1]])
                legend("topright", col=terrain.colors(5), 
                        paste(c(1, defbreaks[1:4]), "-", defbreaks[1:5]), pch=15, 
                        pt.cex=2, cex=0.7, bty="n")
          }})
          
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++ Class: StabilityPCA +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="StabilityPCA",
        representation(eigenvalues = "numeric", measure = "numeric"))
                       
setMethod("show", signature="StabilityPCA", function(object){
          cat("Stability measure: Eigenvalue ratio from Principal Components Analysis, \n",
               "eigenvalue ratio for largest eigenvalue is: ", object@measure, "\n") 
          })
          
          
#+++++++++++ Class: StabilityGLM +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="StabilityGLM",
        representation(coefficients="numeric", deviancevec="numeric",
                       deviancecount = "numeric", weightscheme = "list"))

setMethod("show", signature="StabilityGLM", function(object){
          wl <- object@weightscheme
          if(wl$scheme == "pval") scheme <- "based on p-values"
          else scheme <- "based on ranks"
          decay <- switch(wl$decay, linear="linear weight decay",
                                    quadratic="quadratic weight decay",
                          exponential=paste("exponential weight decay, alpha=", wl$alpha))

          cat("Stability measure: weighted binary logistic regression, \n",
               "weighting: ", scheme, ", ", decay, "\n",
               "deviance criterion is: ", object@deviancecount, "\n")
          })
          
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#+++++++++++ Class: PCAMethodsResult ++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="PCAMethodsResult",
        representation(summary = "numeric", eigenvalues = "numeric", methods = "character"))


setMethod("show", signature="PCAMethodsResult", function(object){
            cat("Aggregation by the first principcal component based on \n ranks for the following methods: \n")
            for(i in 1:length(object@methods)) cat(object@methods[i], "\n")
            cat("variation explained by the first principal component: \n")
            cat(round(object@eigenvalues[1]/sum(object@eigenvalues), 3), "\n")
            })
          

          
            
         

          

          


          



              



                        
                        