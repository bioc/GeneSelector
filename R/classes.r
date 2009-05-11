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
         representation(x = "matrix", y = "factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))

### show-method: can still be improved w.r.t. method output
### -> will be kept for simplicity -- would require 'switch'
                        
setMethod("show", signature="GeneRanking", function(object){
          cat("Ranking by ", object@method, ",\n", sep="")
          cat("number of genes: ", length(object@statistic), ".\n", sep="")
          })

### summary-method

setGeneric("summary", function(object, ...) standardGeneric("summary"))
          
setMethod("summary", signature="GeneRanking", function(object){
          if(is.na(object@pval[1])){ 
           summ <- as.matrix(summary(object@statistic))
           colnames(summ) <- "statistic"
           print(summ)
           }
          else{
                summ <- summary(object@statistic)
                indpp <- match(round(as.numeric(summ),3), round(object@statistic, 3))
                pp <- object@pval[indpp]
                cbind(statistic=summ, p_value = pp)
              }
                                                               }
          )
          
### toplist-method GEÄNDERT 18/11/2008.

setGeneric("toplist", function(object, top=10, show=TRUE,...) standardGeneric("toplist"))

setMethod("toplist", signature(object="GeneRanking"), function(object, top=10, show = TRUE){
          ind <- match(1:top, object@ranking)
          ret <- data.frame(index=ind, statistic=object@statistic[ind], pval=object@pval[ind])
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
              
#+++++++++++ Class: Repeated Ranking +++++++++++++++++++++++++++++++++++++++++++++#
### NOTE: CLASS RENAMED 24/11/2008

setClass(Class="RepeatedRanking",
        representation(original="GeneRanking", rankings="matrix", pvals="matrix",
                       statistics = "matrix", scheme="character"))
                       
setMethod("show", signature="RepeatedRanking", function(object){
          cat(ncol(object@rankings), " rankings with scheme '", object@scheme, "' \n", sep="")
          cat("Method used: ", object@original@method, "\n", sep="")
          })
          
### GEÄNDERT 18/11/2008.
          
setMethod("toplist", signature(object="RepeatedRanking"), function(object, top=10, show=TRUE){
          inds <- apply(object@rankings, 2, function(z) match(1:top, z))  ## ok
          if(top == 1) inds <- t(inds)                                                                                    
          uniqinds <- unique(as.vector(inds))
          emptytab <- numeric(length(uniqinds))
          names(emptytab) <- as.character(uniqinds)
          factinds <- apply(inds, 1, factor, levels=uniqinds)
          genetable <- apply(factinds, 2, function(z){ tab <- table(z) 
                                                       emptytab[names(tab)] <- tab[names(tab)]
                                                       return(emptytab)})
          rown <- rownames(genetable)
          ord <- order(genetable[,1], decreasing=TRUE)
          genetableord <- genetable[ord,,drop = FALSE]
          printgenetable <- data.frame(genetableord)
          rownames(genetableord) <- rown
          colnames(printgenetable) <- paste("Rank ", 1:top, sep="")     

          orig <- match(1:top, object@original@ranking)
          if(object@scheme != "merged (methods)"){
          cat("original dataset: \n \n")  
          print(data.frame(index=orig, statistic=object@original@statistic[orig], pvals=object@original@pval[orig]))
          cat("\n \n")
          toplists <- vector(mode="list", length=ncol(object@rankings))
           
          for(i in 1:ncol(object@rankings)){
          toplists[[i]] <- data.frame(index=inds[,i], statistic=object@statistics[inds[,i],i],
                                      pvals=object@pvals[inds[,i],i])
          if(show){
              cat("Replication", i, ": \n \n")
              print(toplists[[i]])
              cat("\n \n")
            }
           }
          }
          else{
             cat("reference method: \n \n")  
            print(data.frame(index=orig))
          cat("\n \n")
          toplists <- vector(mode="list", length=ncol(object@rankings))
           
          for(i in 1:ncol(object@rankings)){
          toplists[[i]] <- data.frame(index=inds[,i])
          if(show){
              cat("Replication", i, ": \n \n")
              print(toplists[[i]])
              cat("\n \n")
            }
           }
          }  
          #if(show){ 
          cat("In the following table, rownames correspond to gene indices. \n") 
          cat("The columns contain the absolute frequencies for the corresponding ranks \n")
          cat("over all replications. \n") 
          cat("Genes are ordered according to the first column, then to the second, and so on. \n \n")
          print(printgenetable)
          cat("\n \n")
          #}
          return(invisible(toplists))})
          
### GEÄNDERT 24/11/2008.
          
setGeneric("plot")

setMethod("plot", signature("RepeatedRanking", "missing"),
          function(x, frac=1/100, ...){
          if(frac < 0 | frac > 1) stop("'frac' must be between 0 and 1. \n")
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          if(!hasArg(xlab)) ll$xlab <- "Rank in the original dataset"
          if(!hasArg(ylab)) ll$ylab <- "Ranks in perturbed datasets"
          if(!hasArg(cex))  ll$cex <- 0.5
          if(!hasArg(main)) ll$main <- ""
          ngenes <- ceiling(nrow(x@rankings)*frac)
          ll$x <- ll$y <- 1:ngenes 
          ll$type <- "n"
          ll$xlim <- c(1, ngenes)                             
          ll$ylim <- c(0, 2*ngenes)
          do.call("plot", args=ll)
          abline(0,1, lwd=2.5, col="blue")
          ll$type <- "p"
          ngenesind <- match(1:ngenes, x@original@ranking)
          for(i in 1:ncol(x@rankings)){
          ll$y <- x@rankings[ngenesind,i]
          do.call("points", args=ll)               
          }
          }
          )
          
### join-Method replaced by "MergeRanking" 18/11/2008, 24/11/2008
### DOCUMENTATION !

### area of application: same original ranking (-> same method in original ranking),
### but different perturbation schemes (GetRepeatRanking).

setGeneric("MergeRankings", function(RR1, RR2) standardGeneric("MergeRankings"))
                  
setMethod("MergeRankings", signature("RepeatedRanking", "RepeatedRanking"),
           function(RR1, RR2){
           r1 <- RR1@rankings
           r2 <- RR2@rankings
           lr1 <- ncol(r1)
           lr2 <- ncol(r2)                   
          if(!identical(RR1@original, RR2@original))
          stop("Input arguments do not match, different original ranking \n")
          if(nrow(r1) != nrow(r2)) 
          stop("Input arguments do not match, different number of rows in ranking matrix \n")
          if(RR1@original@method != RR2@original@method)
          stop("Input arguments do not match, different ranking methods used \n")
           r12 <- cbind(r1, r2)
           t12 <- cbind(RR1@statistics, RR2@statistics)
           if(any(is.na(RR1@pvals)) || any(is.na(RR2@pvals))) p12 <- matrix(NA, ncol = lr1+lr2,nrow = nrow(r1)) 
           else p12 <- cbind(RR1@pvals, RR2@pvals)
           
           colnames(r12) <- colnames(p12) <- colnames(t12) <- character((lr1+lr2))
           colnames(r12)[1:lr1] <- colnames(p12)[1:lr1] <- colnames(t12)[1:lr1] <- paste(rep(RR1@scheme, lr1), 1:lr1, sep="")
           colnames(r12)[-(1:lr1)] <- colnames(p12)[-(1:lr1)] <- colnames(t12)[-(1:lr1)] <- paste(rep(RR2@scheme, lr2), 1:lr2, sep="")  
           scheme <- "merged (rankings)"
           new("RepeatedRanking", original=RR1@original, rankings=r12, pvals=p12, statistics = t12, scheme=scheme)
          } 
         )
         

### Note 1: 1st method in 'Rlist' will be treated as original.
### Note 2: Aggregated Rankings (over different datasets) are admitted.
### Note 3: Functionality adapted from HeatMapMethods, PCAMethods
### Note 4: HeatMapMethods, PCAmethods will be shortened by referring
### to this method.
### Note 5: Completely different from MergeRankings.
### Note 6: Intention: Aggregation over different test statistics.

setGeneric("MergeMethods", function(Rlist) standardGeneric("MergeMethods"))

setMethod("MergeMethods", signature(Rlist = "list"), function(Rlist) {

if(length(Rlist) < 2) stop("'Rlist' must contain at least two elements \n")
clRlist <- unlist(lapply(Rlist, class))
clind <- (clRlist=="GeneRanking")
if(any(!is.element(clRlist, c("GeneRanking", "AggregatedRanking"))))
stop("All elements of 'Rlist' must be of class 'GeneRanking' or 'AggregatedRanking' \n")

ll <- length(Rlist)
RR <- lapply(Rlist, function(z) slot(z, name="ranking"))
lr <- unlist(lapply(RR, length))
if(length(unique(lr)) != 1) stop("All rankings must have the same length \n")
lmethods <- unlist(lapply(Rlist, function(z) slot(z, name="method")))
if(length(unique(lmethods)) != length(lmethods))
warning("At least one method occurs more than once \n")

nro <- length(RR[[1]])
Rmat <- matrix(unlist(RR), nrow=nro, ncol=ll)
colnames(Rmat) <- lmethods                                                                      
###  NOTE: BY CONVENTION, first entry of the list will constitute the original
###  ranking.
if(class(Rlist[[1]]) != "GeneRanking"){
 original <- new("GeneRanking", x = matrix(0), y = factor(0), statistic = NA,
                  ranking = Rmat[,1], pval = NA, method = slot(Rlist[[1]], "method"),
                  type = "")}
 else original <- Rlist[[1]]


new("RepeatedRanking", original = original, rankings = Rmat[,-1,drop = FALSE],
     pvals = matrix(NA), statistics = matrix(NA), scheme = "merged (methods)")

})


### Convenience function, changed 24/11/2008.
### Note 1: var -> sd. "variance" now rather serves as 'topic name'.
### Note 2: Bootstrap/Jackknife-Variance.
### Note 3: Variance assessement only of 'qualitative nature'. No genuine
### 'estimators' of variance, since replications are far away from i.i.d.
### renamed from 'variance' to 'dispersion', since ranks are fixed, only
### estimators are variant.
### NOTE: 'iqr' does not need a 'center'.

setGeneric("dispersion", function(RR, measure = c("sd", "mad", "iqr"), scheme = c("original", "symmetric", "user"), center = NULL)
            standardGeneric("dispersion"))

setMethod("dispersion", signature("RepeatedRanking"),
           function(RR, measure = c("sd", "mad", "iqr"), scheme = c("original", "symmetric", "user"), center = NULL){
           
            scheme <- match.arg(scheme)
            if(!is.element(scheme, c("original", "symmetric", "user")))
            stop("'scheme' must be one of 'original', 'symmetric' or 'user' \n")
           
            if(scheme == "user" & is.null(center))
            warning("'scheme' = 'user', but 'center' is not specified \n")
           
            if(scheme != "user" & !is.null(center))
            warning("'center' ignored \n")

            measure <- match.arg(measure)
            if(!is.element(measure, c("sd", "mad", "iqr")))
            stop("Invalid 'measure' specified \n")
            
            ###  'iqr' is scheme-independent
            if(measure == "iqr"){
             ranking <- RR@rankings
             iqrs <- apply(ranking, 1, IQR)
             # names(iqrs) <- paste("gene", 1:nrow(ranking))
             # return(iqrs)
            }
            ###
           
           if(scheme == "original"){
            r0 <- RR@original@ranking
            ranking <- RR@rankings
            if(ncol(ranking) < 2) stop("Too few replications in the argument 'RR'; variance cannot be estimated \n")
            if(measure == "sd"){
              sds <- sqrt(1/(ncol(ranking)-1)*rowSums((ranking - r0)^2))
              return(sds)
            }
            else{
               mads <- apply(abs(ranking - r0), 1, median)
               return(mads)
            }
           }
           
           if(scheme == "symmetric"){
              R <- cbind(RR@original@ranking, RR@rankings)
              if(measure == "sd"){
                sds <- apply(R, 1, sd)
                return(sds)
              }
             else{
                ### constant not necessary above, ok.
                mads <- apply(R, 1, mad, constant = 1)
                return(mads)
              }
            }
            
            if(scheme == "user"){

             R <- cbind(RR@original@ranking, RR@rankings)
             if(nrow(R) != length(center))
             stop("Length of 'center' and number of genes disagree \n")
             if(any(center < 0) || any(center > nrow(R)))
             warning("Entries of 'center' negative or larger than the number of genes; verify correctness of 'center' \n")
             
             if(measure == "sd"){
                 sds <- sqrt(1/(ncol(R)-1)*rowSums((R - center)^2))
                 return(sds)
              }
             else{
                mads <- apply(abs(R - center), 1, median)
               return(mads)
              }
            }
          })

          
#+++++++++++ Class: StabilityDistance +++++++++++++++++++++++++++++++++++++++++++++++#
### created 29/11/2008. Replaces StabilityLm.

setClass(Class="StabilityDistance",
        representation(scores ="numeric", noinformation = "numeric", scheme = "character", measure = "character", decay = "character"))
                       
setMethod("show", signature="StabilityDistance", function(object){
          decay <- switch(object@decay, linear = "linear weight decay",
                                    quadratic = "quadratic weight decay",
                                    exponential = "exponential weight decay")
          measure <- switch(object@measure, l1 = "absolute distance",
                                      l2 = "squared distance",
                                      spearman = "spearman's rank correlation",
                                      kendall = "kendall's tau")
          
          cat("Stability measure: ", measure,",\n",
               "scheme: ", object@scheme,",\n",
               "weighting: " , decay ,".\n", sep="")
          })
          
setMethod("summary", signature="StabilityDistance", function(object, display = c("summary", "all"), digits = 3){
           display <- match.arg(display)
           scores <- object@scores
           scheme <- object@scheme
           if(scheme == "original"){
             if(display == "summary"){
               cat("summary of stability scores (with respect to reference data set): \n")
              print(summary(scores, digits = digits))
              cat("expected score in the case of no-information: ", object@noinformation, "\n")
              }
             else{
               cat("stability scores (with respect to reference data set): \n")
               print(round(scores, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation, "\n")
             }
           }
           else{
              if(display == "summary"){
               cat("summary of pairwise stability scores: \n")
               print(summary(scores, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation, "\n")
              }
              else{
                 ### determine number of lists from number of pairs
                 norankings <- 0.5 + sqrt(0.25 + 2*length(scores))
                 scoremat <- matrix(0, nrow = norankings, ncol = norankings)
                 rownames(scoremat) <- colnames(scoremat) <- paste("list", 1:norankings, sep="")
                 diag(scoremat) <- 1
                 k <- 1
                  for(i in 1:norankings){
                    j <- i+1
                    while(j <= norankings){
                    scoremat[i,j] <- scoremat[j,i] <- scores[k]
                    j <- j+1                       
                    k <- k+1
                   }
                  }  
                cat("matrix of pairwise stability scores: \n")
                print(round(scoremat, digits = digits))
                cat("expected score in the case of no-information: ", object@noinformation, "\n")
              }
             }
          })

#+++++++++++ Class: StabilityOverlap +++++++++++++++++++++++++++++++++++++++++++++++#

setClass(Class="StabilityOverlap",
        representation(intersection = "matrix", overlapscore = "matrix",
                       noinformation = "list", scheme = "character", decay = "character"))
        
setMethod("show", signature="StabilityOverlap", function(object){
            decay <- switch(object@decay, linear = "linear weight decay",
                                          quadratic = "quadratic weight decay",
                                          exponential = "exponential weight decay")

            #measure <- switch(object@measure,  intersection = "intersection",
            #                                   union  = "union",
            #                                   overlap_score = "overlap score")

            cat("Stability measure: intersection count and overlap score,\n",
                 "scheme: ", object@scheme, ", \n",
                 "weighting: ", decay, ".\n", sep="")
         })
          
### for generic method definition, s. above
          

setMethod("summary", signature="StabilityOverlap", function(object, measure = c("overlapscore", "intersection"), display = c("summary", "all"), position, digits = 3){
           measure <- match.arg(measure)
           display <- match.arg(display)                                                                                                                                                        
           if(!is.element(measure, c("overlapscore", "intersection")))
           stop("'measure' must be either 'overlapscore' or 'intersection' \n")
           p <- nrow(object@overlapscore)
           if(missing(position)) position <- p
           if(position < 1) stop("'position' must be greater than one \n")
           if(position > p) stop("'position' exceeds the lengths of the lists \n")
           overlapscore <- object@overlapscore[position,]
           intersection <- object@intersection[position,]
           scheme <- object@scheme
           if(scheme == "original"){
             if(display == "summary"){
               if(measure == "overlapscore"){
               cat("summary of overlap scores (with respect to reference data set): \n")                              
               print(summary(overlapscore, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation$overlapscore[position], "\n")
               }
               else{
                  cat("summary of intersection counts (with respect to reference data set): \n")
                  print(summary(intersection, digits = digits))
                  cat("expected score in the case of no-information: ", object@noinformation$intersection[position], "\n")
                }
              }
             else{
               if(measure == "overlapscore"){
               cat("overlap scores (with respect to reference data set): \n")
               print(round(overlapscore, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation$overlapscore[position], "\n")
              }
              else{
               cat("intersection fractions (with respect to reference data set): \n")
               print(round(intersection, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation$intersection[position], "\n")
              }
            }      
           }
           else{
              if(display == "summary"){
               if(measure == "overlapscore"){
               cat("summary of pairwise overlap scores: \n")
               print(summary(overlapscore, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation$overlapscore[position], "\n")
               }
               else{
               cat("summary of pairwise intersection fractions: \n")
               print(summary(intersection, digits = digits))
               cat("expected score in the case of no-information: ", object@noinformation$intersection[position], "\n")
               }
              }
              else{
                 ### determine number of lists from number of pairs.
                 ### NOTE: overlapscore is a vector (s. extraction above).
                 norankings <- as.integer(0.5 + sqrt(0.25 + 2*length(overlapscore)))
                 scoremat <- matrix(0, nrow = norankings, ncol = norankings)
                 rownames(scoremat) <- colnames(scoremat) <- paste("list", 1:norankings, sep="")
                 diag(scoremat) <- 1
                 k <- 1
                 if(measure == "overlapscore"){
                  for(i in 1:norankings){
                    j <- i+1
                    while(j <= norankings){
                    scoremat[i,j] <- scoremat[j,i] <- overlapscore[k]
                    j <- j+1                       
                    k <- k+1
                   }
                  }  
                cat("matrix of pairwise overlap scores: \n")
                print(round(scoremat, digits = digits))
                cat("expected score in the case of no-information: ", object@noinformation$overlapscore[position], "\n")
                }
                else{
                  for(i in 1:norankings){
                    j <- i+1
                   while(j <= norankings){
                    scoremat[i,j] <- scoremat[j,i] <- intersection[k]
                    j <- j+1                     
                    k <- k+1
                    }
                   }  
                  cat("matrix of pairwise intersection fractions: \n")
                  print(round(scoremat, digits = digits))
                  cat("expected score in the case of no-information: ", object@noinformation$intersection[position], "\n")
                }
              }
             }
          })

####
setMethod("plot", signature("StabilityOverlap", "missing"), 
           function(x, frac=1/50, mode = c("mean", "all", "specific"), which = 1, ...){
          dotsCall <- substitute(list(...))
         mode <- match.arg(mode)                                                                              
          ### use recursion
          if(mode == "all"){
            ask <- ((prod(par("mfcol"))) == 1 && dev.interactive())
            opar <- par(ask=ask)

            for(i in 1:ncol(x@overlapscore))
              plot(x, mode = "specific", which = i)
              on.exit(par(opar))                     
          }
          ll <- eval(dotsCall)
          ll$xlab <- "list position"
          ll$ylab <- "overlap"
          ll$lwd <- 1.2
          if(mode == "mean")
          ll$main <- "percentage of  overlap"
          if(mode == "specific")
          ll$main <- paste("percentage of overlap, which", which, sep = "=")
          ll$type <- "l"
          ll$col <- "grey"
          ngenes <- ceiling(nrow(x@overlapscore)*frac)
          ll$x <- 1:ngenes
          if(mode == "mean")
            ll$y <- apply(x@intersection[1:ngenes, , drop=FALSE], 1, mean)
          if(mode == "specific")
            ll$y <- x@intersection[1:ngenes, which]
          ll$ylim <- c(0, 1)  ### changed
          layout(mat=as.matrix(c(1,2)), heights=c(1,1))
          do.call("plot", args=ll)
          abline(h = 1, lwd=2.5)
          ll$col <- "red"
          ## ll$lty <- "dashed"
          ll$lwd <- 1.2
          ll$type <- "l"
          if(mode == "mean")
            ll$y <- apply(x@overlapscore[1:ngenes, , drop=FALSE], 1, mean)
          if(mode == "specific")
            ll$y <- x@overlapscore[1:ngenes, which]
          ll$ylim <- c(0,1)
          ll$ylab <- "score"
          if(mode == "mean")
            ll$main <- "average overlap score"
          else
            ll$main <- paste("overlap score, which", which, sep = "=")
          do.call("plot", args=ll)
          abline(h = 1, lwd = 2.5)                                                                             
          layout(mat=as.matrix(1))                                                                                
           }
          )
          
          
####

setClass(Class="StabilityUnion",
        representation(union = "numeric", unionscore = "numeric", noinformation = "list", decay = "character"))

setMethod("show", signature="StabilityUnion", function(object){
            cat("Stability measure: union count, \n", sep="")
            cat("weighting: " , object@decay , ". \n", sep ="")
         })
         
### plot-method


setMethod("plot", signature("StabilityUnion", "missing"),
           function(x, frac=1/50, ...){
          dotsCall <- substitute(list(...))
          ll <- eval(dotsCall)
          ll$xlab <- "list position"
          ll$ylab <- "union"
          ll$lwd <- 1.2
          ll$main <- paste("union")
          ll$type <- "l"
          ll$col <- "grey"
          ngenes <- ceiling(length(x@unionscore)*frac)
          ll$x <- 1:ngenes
          ll$y <- x@union[1:ngenes]
          ll$ylim <- c(0, 1)
          layout(mat=as.matrix(c(1,2)), heights=c(1,1))
          do.call("plot", args=ll)
          abline(h = 1, lwd = 2.5)                             
          if(length(x@noinformation) != 0){
            ll$lty <- "dashed"
            ll$y <- x@noinformation$union[1:ngenes]
            do.call("lines", args = ll)
          }
          #abline(0,1, lwd=2.5)
          ll$col <- "red"
          ll$lwd <- 1.2
          ll$type <- "l"
          ll$y <- x@unionscore[1:ngenes]
          ll$ylim <- c(0,1)
          ll$lty <- 1                             
          ll$ylab <- "score"
          ll$main <- "union score"
          do.call("plot", args=ll)
          abline(h = 1, lwd = 2.5)                             
          if(length(x@noinformation) != 0){
            ll$y <- x@noinformation$unionscore[1:ngenes]
            ll$lty <- "dashed"                               
            do.call("lines", args = ll)
          }

          layout(mat=as.matrix(1))
           }
          )


          
         
#+++++++++++ Class: AggregatedRanking +++++++++++++++++++++++++++++++++++++++++#

### major revision 18/11/2008.
### type: "which aggregation method" (simple, MC, PCA, ...)
### measure: "Was war Grundlage der Aggregation (entsprich Psi_j in den Definitionen)"
### method: Welche Methode wurde aggregiert (nicht immer sinnvoll).

setClass(Class="AggregatedRanking",
        representation(ranking = "numeric", type= "character", measure = "character", method = "character"))
        
           
setMethod("show", signature="AggregatedRanking", function(object){
            cat(object@type, " aggregation for ", length(object@ranking), "genes, \n", sep="")
            cat("aggregation measure: ",  object@measure, ".\n", sep="")
            })

setMethod("toplist", signature(object="AggregatedRanking"), function(object, top=10, show = TRUE){
          ind <- match(1:top, object@ranking)
          ret <- data.frame(index=ind)
          rownames(ret) <- NULL
          if(show) print(ret)
          invisible(ret)
          })


#+++++++++++ Class: GeneSelectorOutput +++++++++++++++++++++++++++++++++++++++++#
### Result after running GeneSelector.
### class renamed 17/12/2008.

setClass(Class="GeneSelectorOutput",
        representation(final = "numeric", rankings = "matrix", inout = "matrix", 
                       selected = "numeric", adjpval="vector",  maxrank = "numeric", 
                       statistics = "character"))
        
           
setMethod("show", signature="GeneSelectorOutput", function(object){
            cat("GeneSelector run with gene rankings from the following statistics: \n")
            for(i in 1:length(object@statistics)) cat(object@statistics[i], "\n")
            cat("Number of genes below threshold rank ", 
                 object@maxrank, " in all statistics:", 
                 sum(object@selected), "\n", sep="")})
                 
setMethod("toplist", signature=(object="GeneSelectorOutput"), function(object, top=10, show = TRUE){
          ind <- match(1:top, object@final)
          ret <- data.frame(index=ind, pvals=object@adjpval[ind])
           rownames(ret) <- NULL
          if(show) print(ret)
          invisible(ret)
          })
          
setGeneric("SelectedGenes", function(object) standardGeneric("SelectedGenes"))

setMethod("SelectedGenes", signature(object="GeneSelectorOutput"), function(object){
          toplist(object, sum(object@selected)) })
          

setMethod("plot", signature("GeneSelectorOutput", "missing"),
          function(x, which){
          if(length(which) != 1) stop("Length of 'which' must be one. \n")
          stats <- rev(x@statistics)
          nostats <- length(stats)
          inout <- rev(x@inout[which,])
          RR <- rev(x@rankings[which,])
          p <- nrow(x@rankings)
          p0 <- c(0.01, 0.05, 0.1, 0.25, 1)
          defbreaks <- ceiling(p0*p)
          plot(0:1, 0:1, axes=FALSE, type="n", 
              main=paste("GeneInfoScreen for gene ", substitute(which)),
              xlab="", ylab="")
          partition <- seq(from=0, to=0.9, length=nostats+1)
          cols <- terrain.colors(5)
          for(i in 1:nostats){
          rect(0, partition[i], 0.25, partition[i+1])
          points(0.125, (partition[i] + partition[i+1])/2, cex=1.7, pch=inout[i])
          text(0.15, 1, "selected ?", cex=1)
          text(0.5, 1,  "criterion", cex=1)    ### statistic -> criterion.
          text(0.8, 1, "rank", cex=1)
          text(0.5, (partition[i] + partition[i+1])/2, stats[i], cex=1.5)
          text(0.8, (partition[i] + partition[i+1])/2, round(RR[i]), cex=1.5,
               col = cols[which(defbreaks >= RR[i])[1]])
                legend("topright", col=terrain.colors(5), 
                        paste(c(1, defbreaks[1:4]), "-", defbreaks[1:5]), pch=15, 
                        pt.cex=2, cex=0.7, bty="n")
          }})
          
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


          



              



                        
                        
