### filename: RecoveryScore.r
### Title:
### Compute Recover Score of Pavlidis et al. (2003)
###
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 10.9.2007
### date(s) of updates: 11.9.2007
#
### Brief description:
#
#  Computes Recovery Scores for objects of class
#  of an object of class RepeatRanking that
#  must contain pvalues for each replication.
#  Pvalues can be adjusted using the function 'AdjustPvalues'.
#
#
#

#
### Further comments and notes:
#
#
#
#
#
#
###**************************************************************************###

RecoveryScore <- function(RR, method=c("raw", "BH", "qvalue", "Bonferroni",
                              "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"),
                                maxpval=0.05){
 if(class(RR) != "RepeatRanking") 
 stop("'RR' must be an object of class 'RepeatedRanking' \n")                               
 oP <- RR@original@pval
 oR <- RR@original@ranking
 P <- RR@pvals
 R <- RR@rankings
 if(any(is.na(oP)) | any(is.na(P)))
 stop("For the computation of the Recovery Score complete sets
       of p-values must be given \n")
 if(method == "raw"){ 
 oadjP <- oP
 adjP <- P
 }
 else{ 
 oadjP <- AdjustPvalues(oP, method=method)
 adjP <- apply(P, 2, function(z) AdjustPvalues(z, method=method))
 }
 ocountsignif <- length(oadjP[oadjP<=maxpval])
 osel <- oR[1:ocountsignif]
 countsignif <- apply(adjP, 2, function(z) length(z[z<=maxpval]))
 recoveryscores <- double(ncol(R))
 for(i in 1:ncol(R)){
 recoveryscores[i] <- length(intersect(R[1:countsignif[i],i], osel))/countsignif[i]
 if(!is.finite(recoveryscores[i])) {
  recoveryscores[i] <- 0
  warning("No significant genes found ; recovery score set to 0. \n")
  }
 }
 names(recoveryscores) <- as.character(1:ncol(R))
 return(recoveryscores)
 }
 

                              
                              
                              
                              
