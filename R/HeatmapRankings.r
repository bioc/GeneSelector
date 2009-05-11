### filename: HeatmapMethods.r
### Title: Create a heatmap for various gene ranking statistics.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 28.8.2007
### date(s) of updates: 7.9.2007, 12.9.2007
### major revision: 17.12.2008 
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

setGeneric("HeatmapRankings", function(RR, ind=1:100) standardGeneric("HeatmapRankings"))

setMethod("HeatmapRankings", signature(RR = "RepeatedRanking"), function(RR, ind=1:100) {

Rmat <- cbind(RR@original@ranking, RR@rankings)
rownames(Rmat) <- rownames(RR@original@ranking)
colnames(Rmat) <- c(RR@original@method, colnames(RR@rankings))
p0 <- c(0.01, 0.05, 0.1, 0.25, 1)
defbreaks <- ceiling(p0*nrow(Rmat))
hm <- heatmap(Rmat[ind, , drop = FALSE], hclustfun = hclust, labCol=colnames(Rmat), scale="none",
              labRow = character(), col = terrain.colors(20), 
              mar=c(8,8))
legend("topright", col=terrain.colors(5), 
        legend=paste(c(1, defbreaks[1:4]), "-", defbreaks[1:5]), pch=15, 
        pt.cex=2, cex=0.7, bty="n")
return(invisible(hm))

})







 


















            
            
            
            
