### filename: PCAMethods.r
### Aggregate over different statistics using PCA.
###**************************************************************************###

setGeneric("PCAMethods", function(Rlist) standardGeneric("PCAMethods"))

setMethod("PCAMethods", signature(Rlist="list"), function(Rlist) {

if(length(Rlist) < 2) stop("'Rlist' must contain at least two elements \n")
clRlist <- unlist(lapply(Rlist, class))
clind <- (clRlist=="GeneRanking")
if(any(!is.element(clRlist, c("GeneRanking", "AggregatedRanking"))))
stop("All Elements of 'Rlist' must be of class 'GeneRanking' or 'AggregatedRanking' \n")

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
omat <- matrix(nrow=nro, ncol=ll)
if(any(clRlist == "GeneRanking"))
omat[,clind] <- matrix(apply(xmat[,clind,drop=FALSE], 2, function(z) match(1:nro, z)), ncol=ll)
if(any(clRlist == "AggregatedRanking"))
omat[,!clind] <- xmat[,!clind]
pcares <- prcomp(omat, retx = TRUE, center = TRUE)
rot <- pcares$x[,1]
ordrot <- order(rot)
agg <- match(1:lr[1], ordrot)
eigenval <- (pcares$sdev)^2
new("PCAMethodsResult", summary = agg, eigenvalues = eigenval, methods = lmethods)

})