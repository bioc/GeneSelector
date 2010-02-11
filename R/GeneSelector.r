### filename: GeneSelector.r
### Title: Create a set of 'surviving genes' using
### a user defined ranking of ranking methods
### and a threshold value either user-defined or
### estimated via adjusted p-values/number of differentially
### expressed genes.
###
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 7.9.2007
### date(s) of updates: 10.9.2007, 11.9.2007
#
### Brief description:
#
#  Returns a final ranking/list of 'surviving' genes,
#  based upon single rankings or aggregated repeated rankings.
#
#

#
### Further comments and notes:
#
#
#
###**************************************************************************###

setGeneric("GeneSelector", function(Rlist, ind=NULL, indstatistic=1:length(Rlist),
threshold=c("user", "BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS",
        "SidakSD", "BY"), maxrank=NULL, maxpval=0.05)
standardGeneric("GeneSelector"))

setMethod("GeneSelector", signature(Rlist="list"), function(Rlist, ind=NULL, indstatistic=1:length(Rlist),
                     threshold=c("user", "BH", "qvalue", "Bonferroni", "Holm",
                     "Hochberg", "SidakSS", "SidakSD", "BY"), maxrank=NULL, maxpval=0.05) {

threshold <- match.arg(threshold)
if(!is.element(threshold, eval(formals(GeneSelector)$threshold)))
stop("Invalid threshold method defined. \n")
if(threshold == "user" & is.null(maxrank))
stop("'maxrank' must be specified if 'threshold' is 'user'")
if(threshold != "user" & (maxpval <= 0 | maxpval > 1))
stop("'maxpval' must between 0 and 1")
clRlist <- unlist(lapply(Rlist, class))
clind <- (clRlist == "GeneRanking")
if(any(!is.element(clRlist, c("GeneRanking", "AggregatedRanking"))))
stop("All Elements of 'Rlist' must be of class 'GeneRanking'
      or 'AggregatedRanking' \n")
ll <- length(Rlist)
RR <- lapply(Rlist, function(z)  slot(z, name = "ranking"))
lr <- unlist(lapply(RR, length))
if(length(unique(lr)) != 1) stop("All Rankings must have the same length \n")
lmethods <- unlist(lapply(Rlist, function(z) slot(z, name="method")))
if(length(unique(lmethods)) != length(lmethods))
warning("One or more methods occur more than once \n")
nro <- length(RR[[1]])
if(is.null(maxrank)){
pvallist  <- lapply(Rlist, function(z) slot(z, name="pval"))
for(k in 1:length(pvallist)){
 if(any(is.na(pvallist[[k]]))) next
 else { pval <- pvallist[[k]] ; break }
 }
 if(inherits(try(pval <- pval), "try-error"))
  stop("None of the provided rankings contains pvalues required for
        threshold finding \n")
  else{
       pval <- AdjustPvalues(pval, threshold)
       maxrank <- length(pval[pval<=maxpval])
      }
 }
else pval <- rep(NA, nro)
xmat <- matrix(unlist(RR), nrow=nro, ncol=ll)
if(is.null(ind)) ind <- 1:nro
colnames(xmat) <- lmethods
indstatistic <- as.integer(indstatistic)
if(any(indstatistic <= 0 | indstatistic > ll))
stop("Invalid specification of 'indstatistic' \n")
xmatord <- xmat[,indstatistic,drop=FALSE]
xmatbool <- 1*(xmatord <= maxrank)
inout <- apply(xmatbool, 2, function(z) ifelse(z == 1, "+", "-"))
survind <- (rowSums(xmatbool) == length(indstatistic))
arglistord <- data.frame(cbind(1*(!xmatbool), xmatord)) 
ranking <- do.call("order", arglistord)
#absdist <- rowSums(xmat) - ll
#maxdist <- (nro-1)*ll
#reldist <- absdist/maxdist 
new("GeneSelectorOutput", final = match(1:nro, ranking), rankings = xmatord,  inout = inout,  
    selected = as.numeric(survind), adjpval = pval, maxrank = maxrank,
    statistics = lmethods[indstatistic])
})

