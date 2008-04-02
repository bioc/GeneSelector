### filename: HeatmapMethods.r
### Title: Create a heatmap for various gene ranking statistics.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.8.2007
### date(s) of updates: 7.9.2007, 12.9.2007
#
### Brief description:
#
#  Generates a heatmap with
#  columns = ranking methods
#  rows = genes
#  from hierarchical (complete linkage) clustering for both
#  genes and methods.
#

#
### Further comments and notes:
#
#   Explorative method.
#
###**************************************************************************###

setGeneric("HeatmapMethods", function(Rlist, ind=1:100) standardGeneric("HeatmapMethods"))

setMethod("HeatmapMethods", signature(Rlist="list"), function(Rlist, ind=1:100) {

if(length(Rlist) < 2) stop("'Rlist' must contain at least two elements \n")
clRlist <- unlist(lapply(Rlist, class))
clind <- (clRlist=="GeneRanking")
if(any(!is.element(clRlist, c("GeneRanking", "AggregatedRanking"))))
stop("All Elements of 'Rlist' must be of class 'GeneRanking' 
      or 'AggregatedRanking' \n")

ll <- length(Rlist)
RR <- lapply(Rlist, function(z){
             namz <- slotNames(z)   
             if(is.element("ranking", namz)) slot(z, name="ranking")
             else slot(z, name = "summary")})
lr <- unlist(lapply(RR, length))
if(length(unique(lr)) != 1) stop("All Rankings must have the same length \n")
lmethods <- unlist(lapply(Rlist, function(z) slot(z, name="method")))
if(length(unique(lmethods)) != length(lmethods))
warning("One or more methods occur more than once \n")

nro <- length(RR[[1]])
xmat <- matrix(unlist(RR), nrow=nro, ncol=ll)
### change to ranks:
if(is.null(ind)) ind <- 1:nro
omat <- matrix(nrow=length(ind), ncol=ll)
if(any(clRlist == "GeneRanking"))
omat[,clind] <- matrix(apply(xmat[,clind,drop=FALSE], 2, function(z) match(ind, z)), ncol=sum(clind))
if(any(clRlist == "AggregatedRanking"))
omat[,!clind] <- xmat[ind,!clind]
rownames(omat) <- rownames(RR[[1]])
colnames(omat) <- lmethods
p0 <- c(0.01, 0.05, 0.1, 0.25, 1)
defbreaks <- ceiling(p0*nro)
hm <- heatmap(omat, hclustfun = hclust, labCol=lmethods, scale="none",
              labRow = character(), col = terrain.colors(20), 
              mar=c(8,8))
legend("topright", col=terrain.colors(5), 
        legend=paste(c(1, defbreaks[1:4]), "-", defbreaks[1:5]), pch=15, 
        pt.cex=2, cex=0.7, bty="n")
return(invisible(hm))

})







 


















            
            
            
            