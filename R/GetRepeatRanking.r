### filename: GetRepeatRanking.r
### Title: Generate repeat rankings from perturbed datasets.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 17.8.2007
### date(s) of updates: 30.8.2007
#
### Brief description:
#
#   Generates repeated rankings in three possible manners:
#   - Subsampling (with constraints y/n)
#   - switching class labels
#   - Bootstrapping (with constraints y/n)
#   - Adding Gaussian Noise
#

#
### Further comments and notes:
#
#   s. also GenerateFoldMatrix.r, GenerateBootMatrix.r
#
###**************************************************************************###

################################################################################

setClass(Class="BootMatrix",
        representation(bootmatrix="matrix", replicates="integer",
                       type="character", maxties="vector", minclassize="integer", 
                        balancedclass="logical", balancedsample="logical"))

setClass(Class="FoldMatrix",
        representation(foldmatrix="matrix", k="integer", replicates="integer",
                       type="character", minclassize="integer", balanced="logical"))


setClass(Class="GeneRanking",
        representation(x="matrix", y="factor", statistic="numeric", ranking="numeric",
                        pval="vector", type="character", method="character"))

################################################################################

setGeneric("GetRepeatRanking", function(R, P, scheme=c("Subsampling", "Labelexchange"),
                                        iter=10, varlist=list(genewise=FALSE, factor=1/5), ...)
             standardGeneric("GetRepeatRanking"))
             
### Subsampling:

setMethod("GetRepeatRanking", signature(R="GeneRanking", P="FoldMatrix", iter="missing", varlist="missing"),
            function(R, P, scheme=c("Subsampling", "Labelexchange"), ...){
            scheme <- match.arg(scheme)
            if(!is.element(scheme, c("Subsampling", "Labelexchange")))
            stop("'scheme' must be  either 'Subsampling' or 'Labelexchange'")
            x <- R@x
            y <- R@y
            Pm <- P@foldmatrix
            type <- R@type
            iter <- ncol(Pm)
            rankm <- pvalm <- statisticm <- matrix(nrow=nrow(x), ncol=iter)
            rankfun <- switch(R@method, ordinaryT=RankingTstat,
                                        WelchT=RankingWelchT,
                                        BaldiLongT=RankingBaldiLong,
                                        Bstat=RankingBstat,
                                        Ebam=RankingEbam,
                                        Foldchange=RankingFC,
                                        FoxDimmicT=RankingFoxDimmic,
                                        Gapstatistic=RankingGap,
                                        Limma=RankingLimma,
                                        Permutation=RankingPermutation,
                                        Sam=RankingSam,
                                        ShrinkageT=RankingShrinkageT,
                                        SoftthresholdT=RankingSoftthresholdT,
                                        WilcEbam=RankingWilcEbam,
                                        Wilcoxon=RankingWilcoxon)
            if(scheme == "Subsampling"){
            for(i in 1:iter){
             currx <- x[,Pm[,i]]
             curry <- y[Pm[,i]]
             repet <- rankfun(currx, curry, type, ...)
             rankm[,i] <-  repet@ranking
             pvalm[,i] <- repet@pval
             statisticm[,i] <- repet@statistic
            }
            }
            if(scheme == "Labelexchange"){
            ly <- levels(y)
            nly <- nlevels(y)
            if(nly != 2) stop("scheme 'Labelexchange' not allowed if y has only one level \n")
            for(i in 1:iter){
             curry <- y
             curry[!Pm[,i]] <- ifelse(y[!Pm[,i]] == ly[1], ly[2], ly[1])
             repet <- rankfun(x, curry, type, ...)
             rankm[,i] <-  repet@ranking
             pvalm[,i] <- repet@pval
             statisticm[,i] <- repet@statistic
            }
            }
            
          colnames(rankm) <- colnames(pvalm) <- colnames(statisticm) <- paste("iteration", 1:iter)
          rownames(rankm) <- rownames(pvalm) <- rownames(statisticm) <- paste("top gene", 1:nrow(x)) 
            
            new("RepeatRanking", original=R, rankings=rankm, pvals=pvalm,
                statistics=statisticm, scheme=scheme)
            }
            )
            
### Bootstrap

setMethod("GetRepeatRanking", signature(R="GeneRanking", P="BootMatrix", scheme = "missing",
                                        iter = "missing", varlist = "missing"),
           function(R, P,...){
            x <- R@x
            y <- R@y
            Pm <- P@bootmatrix
            type <- R@type
            iter <- ncol(Pm)
            rankm <- pvalm <- statisticm <- matrix(nrow=nrow(x), ncol=iter)
            rankfun <- switch(R@method, ordinaryT=RankingTstat,
                                        WelchT=RankingWelchT,
                                        BaldiLongT=RankingBaldiLong,
                                        Bstat=RankingBstat,
                                        Ebam=RankingEbam,
                                        Foldchange=RankingFC,
                                        FoxDimmicT=RankingFoxDimmic,
                                        Gapstatistic=RankingGap,
                                        Limma=RankingLimma,
                                        Permutation=RankingPermutation,
                                        Sam=RankingSam,
                                        ShrinkageT=RankingShrinkageT,
                                        SoftthresholdT=RankingSoftthresholdT,
                                        WilcEbam=RankingWilcEbam,
                                        Wilcoxon=RankingWilcoxon)
            for(i in 1:iter){
             currx <- x[,Pm[,i]]
             curry <- y[Pm[,i]]
             repet <- rankfun(currx, curry, type, ...)
             rankm[,i] <-  repet@ranking
             pvalm[,i] <- repet@pval
             statisticm[,i] <- repet@statistic
            }
            
          colnames(rankm) <- colnames(pvalm) <- colnames(statisticm) <- paste("iteration", 1:iter)
          rownames(rankm) <- rownames(pvalm) <- rownames(statisticm) <- paste("top gene", 1:nrow(x)) 

            new("RepeatRanking", original=R, rankings=rankm, pvals=pvalm,
                statistics=statisticm, scheme="Bootstrap")
           })

### Adding noise

setMethod("GetRepeatRanking", signature(R="GeneRanking", P="missing", scheme = "missing"),
           function(R, iter=10, varlist=list(genewise=FALSE, factor=1/5), ...){
           genewise <- varlist$genewise
           if(is.null(genewise)) genewise <- FALSE
           factor <- varlist$factor
           if(is.null(factor)) factor <- 1/5
            x <- R@x
            y <- R@y
            ly <- length(y)
            type <- R@type
            rankm <- pvalm <- statisticm <- matrix(nrow=nrow(x), ncol=iter)
            rankfun <- switch(R@method, ordinaryT=RankingTstat,
                                        WelchT=RankingWelchT,
                                        BaldiLongT=RankingBaldiLong,
                                        Bstat=RankingBstat,
                                        Ebam=RankingEbam,
                                        Foldchange=RankingFC,
                                        FoxDimmicT=RankingFoxDimmic,
                                        Gapstatistic=RankingGap,
                                        Limma=RankingLimma,
                                        Permutation=RankingPermutation,
                                        Sam=RankingSam,
                                        ShrinkageT=RankingShrinkageT,
                                        SoftthresholdT=RankingSoftthresholdT,
                                        WilcEbam=RankingWilcEbam,
                                        Wilcoxon=RankingWilcoxon)
            if(!genewise){
              sigma <- factor*sd(x)
              for(i in 1:iter){
              jittering <- matrix(rnorm(prod(dim(x)), mean=0, sd=sigma),
                                  nrow=nrow(x), ncol(x))
              currx <- x+jittering
              repet <- rankfun(currx, y, type, ...)
              rankm[,i] <-  repet@ranking
              pvalm[,i] <- repet@pval
              statisticm[,i] <- repet@statistic
              }
            }
            else{
              sigmavec <- apply(x, 1, sd)*factor
              for(i in 1:iter){
              jittering <- t(sapply(sigmavec, function(z) rnorm(ly, mean=0, sd=z)))
              currx <- x+jittering
              repet <- rankfun(currx, y, type, ...)
              rankm[,i] <-  repet@ranking
              pvalm[,i] <- repet@pval
              statisticm[,i] <- repet@statistic
              }
            }
          
          colnames(rankm) <- colnames(pvalm) <- colnames(statisticm) <- paste("iteration", 1:iter)
          rownames(rankm) <- rownames(pvalm) <- rownames(statisticm) <- paste("top gene", 1:nrow(x)) 
             
          new("RepeatRanking", original=R, rankings=rankm, pvals=pvalm,
                statistics=statisticm, scheme="Jittering")
        }
       )


















