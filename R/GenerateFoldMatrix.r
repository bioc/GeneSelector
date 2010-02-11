### filename: GenerateFoldMatrix.r
### Title: Genrate boolean matrices that define subsampling.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.8.2007
### date(s) of updates:
#
### Brief description:
#
#   Generates a matrix that is used for subsampling and
#   Repeated Ranking of Genes. Subsampling is done
#   with constraints (unique replications, minclassize etc.)
#   Up to now, algorithm is not determistic.
#
#
### Further comments and notes:
                                        #
#   Very difficult and unperfect function.
#   Checks needed.
#   Can take a lot of time for big samplesize/k.
#   See also BootMatrix.r.
#
###**************************************************************************###

### generic function

setGeneric("GenerateFoldMatrix", function(x, y, k=1, replicates = ifelse(k==1, length(y), 10), type=c("unpaired", "paired", "onesample"),
            minclassize=2, balanced=FALSE, control) standardGeneric("GenerateFoldMatrix"))

### signature: x=matrix, y=numeric.

setMethod("GenerateFoldMatrix", signature(x = "missing", y="numeric"),
          function(y, k=1, replicates=ifelse(k==1, length(y), 10), type=c("unpaired", "paired", "onesample"),
         minclassize=2, balanced=FALSE, control){
type <- match.arg(type)
k <- as.integer(k)
ly <- length(y)
if(missing(control)) control <- samplingcontrol(candreplicates=replicates*3)
if (type == "paired" & (k<1 | k > ly/2-1)) stop("Invalid k (<1 or > (length of pairs - 1). \n")
else if(k < 1 | k > ly-2) stop("Invalid k (<1 or > (length of samples - 2)). \n")
  if( !is.element(type, eval(formals(GenerateFoldMatrix)$type)))
  stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
  taby <- table(y)
  if(length(taby) != 2 & type != "onesample")
  stop("y has not exactly two levels, but 'type' is not 'onesample'. \n")
  if(length(taby) != 1 & type == "onesample")
  warning("y has not exactly one level, but 'type' is 'onesample'. \n")
  if(type == "paired"){
  if(length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
           stop("Incorrect coding for type='paired'. \n")
  }
                                                 if(k == 1){
    if(is.element(type, c("unpaired", "onesample"))){
      foldm <- diag(1, nrow=ly)
      mode(foldm) <- "logical"
      foldm <- !foldm
      if(replicates < ly){ 
        foldm <- foldm[,sample(ncol(foldm), replicates),drop=FALSE]
      }
      else replicates <- ly
      type <- type
      minclassize <- min(taby) - 1
      balanced <- FALSE
      }
    else{
      foldm <- matrix(nrow=ly, ncol=ly/2, data=1)
      for(i in 1:(ly/2)) foldm[i,i] <- foldm[ly/2+i, i] <- 0
      mode(foldm) <- "logical"
      minclassize <- ly/2 - 1
      balanced <- FALSE
      if(replicates < ly/2){ 
        foldm <- foldm[,sample(ncol(foldm), replicates),drop=FALSE]
      }
      else replicates <- ly/2
      }
    }
    else{
      replicates <- as.integer(replicates)
      if(replicates < 1) stop("replicates must not be smaller than 1 ! \n")
      if(type == "unpaired"){
        if(balanced) minclassize <- floor((ly - k)/2)
        if(minclassize < 1){
         minclassize <- 1
         warning("Invalid minclassize for type='unpaired'. Reset to 1 \n")
        }
        if((min(table(y))-k) >= minclassize){
         checkfun <- function(z){
         row <- rep(NA,ly)
         row[z] <- y[z]
         return(row)
         }
        }
        else{
         checkfun <- function(z){
         row <- rep(NA,ly)
         row[z] <- y[z]
         tab <- table(row)
         y1 <- tab[names(taby)[1]]
         if(is.na(y1)) y1 <- 0                         
         y2 <- tab[names(taby)[2]]
         if(is.na(y2)) y2 <- 0                        
         if((taby[1]-y1) >= minclassize && (taby[2]-y2) >= minclassize) return(row)
         else return(rep(NA,ly))}
         }
        if(choose(ly,k) > 10000){
        iter <- 0
        UNIQUE <- FALSE
        candreplicates <- with(control, candreplicates)
        maxiter <- with(control, maxiter)
        while(iter < maxiter && !UNIQUE){
        foldm <- replicate(candreplicates, sample(ly, k))
        iter <- iter+1
        foldm <- apply(foldm, 2, checkfun)
        indfoldm <- apply(foldm, 2, function(z) any(!is.na(z)))
        if(length(indfoldm) < replicates){
          if(iter < maxiter){
          cat("Not enough appropriate replications found. Repeating... \n")
          next
          }
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity. \n")
          uniqfoldm <- unique(foldm[,indfoldm], MARGIN=2)
          replicates <- ncol(uniqfoldm)
          }
          }
        else{
        foldm <- foldm[,indfoldm]
        uniqfoldm <- unique(foldm, MARGIN=2)
        UNIQUE <- (ncol(foldm) >= replicates)
        if(!UNIQUE){
            if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
            else{
            warning("Desired Number of replications could not be generated. \n Check all arguments for validity. \n")
            replicates <- ncol(foldm)
            }
           }
          }
         }
         foldm <- is.na(uniqfoldm)
        }
        else{
        foldm <- combn(ly,k)
        foldm <- apply(foldm, 2, checkfun)
        indfoldm <- apply(foldm, 2, function(z) any(!is.na(z)))
        foldm <- foldm[,indfoldm]
        nsc <- ncol(foldm)
        if(replicates > nsc){
        replicates <- nsc
        warning("Desired Number of replications could not be generated. \n Either reduce the number of replicates or change 'minclassize'/'balance' arguments. \n")
        }
        foldm <- is.na(foldm[,,drop=FALSE])
        }
        }
        if(type == "onesample"){
          if(choose(ly,k) > 10000){
        UNIQUE <- FALSE
        iter <- 0
        candreplicates <- with(control, candreplicates)
        maxiter <- with(control, maxiter)
        while(iter < maxiter && !UNIQUE){
        foldm <- replicate(candreplicates, sample(ly, k))
        iter <- iter+1
        foldm <- apply(foldm, 2, function(z){ row <- rep(NA,ly); row[z] <- 1; return(row)})
        uniqfoldm <- unique(foldm, MARGIN=2)
        UNIQUE <- (ncol(uniqfoldm) >= replicates)
        if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications found. Repeating... \n")
          else{
           warning("Desired Number of replications could not be generated. \n Check all arguments for validity. \n")
           replicates <- ncol(foldm)
          }
         }
        }
        foldm <- is.na(uniqfoldm)
       }
        else{
        foldm <- combn(ly,k)
        foldm <- apply(foldm, 2, function(z){ row <- rep(NA,ly); row[z] <- 1; return(row)})
        nsc <- ncol(foldm)
        if(replicates > nsc){
        replicates <- nsc
        warning("Desired Number of replications could not be generated. \n Either reduce the number of replicates or change 'minclassize'/'balance' arguments. \n")
        }
        foldm <- is.na(foldm[,,drop=FALSE])
        }
        balanced <- FALSE
        minclassize <- 0
        }
        if(type == "paired"){
          if(choose(ly/2,k) > 10000){
        UNIQUE <- FALSE
        iter <- 0
        candreplicates <- with(control, candreplicates)
        maxiter <- with(control, maxiter)
        while(iter < maxiter && !UNIQUE){
        foldm <- replicate(candreplicates, sample(ly/2, k))
        iter <- iter+1
        foldm <- apply(foldm, 2, function(z){ row <- rep(NA,ly/2); row[z] <- 1; return(row)})
        uniqfoldm <- unique(foldm, MARGIN=2)
        UNIQUE <- (ncol(foldm) >= ncol(uniqfoldm))
        if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications found. Repeating... \n")
          else{
           warning("Desired Number of replications could not be generated. \n Check all arguments for validity. \n")
           replicates <- ncol(foldm)
          }
         }
        }
        foldm <- is.na(uniqfoldm)
       }
        else{
        foldm <- combn(ly/2,k)
        foldm <- apply(foldm, 2, function(z){ row <- rep(NA,ly/2); row[z] <- 1; return(row)})
        nsc <- ncol(foldm)
        if(replicates > nsc){
        replicates <- nsc
        warning("Desired Number of replications could not be generated. \n Check all arguments for validity. \n")
        }
        foldm <- is.na(foldm[,,drop=FALSE])
        }
        foldm <- rbind(foldm, foldm)
        balanced <- TRUE
        minclassize <- ly/2
      }
     if(ncol(foldm)>replicates) foldm <- foldm[,sample(ncol(foldm), replicates),drop=FALSE]
   }
   new("FoldMatrix", foldmatrix=foldm, k=k, replicates=as.integer(replicates),
                       type=type, balanced=balanced, minclassize=as.integer(minclassize))
  }
)

### signature: x=matrix, y=factor.

setMethod("GenerateFoldMatrix", signature(x = "missing", y="factor"),
function(y, k=1, replicates=ifelse(k==1, length(y), 10), type=c("unpaired", "paired", "onesample"),
         minclassize=2, balanced=FALSE, control){
         GenerateFoldMatrix(as.numeric(y),k,replicates,type,minclassize, balanced,
                           control)})
                           
### signature: x=ExpressionSet, y=character

setMethod("GenerateFoldMatrix", signature(x="ExpressionSet", y="character"),
function(x, y, k=1, replicates=ifelse(k==1, ncol(exprs(x)), 10), type=c("unpaired", "paired", "onesample"),
         minclassize=2, balanced=FALSE, control){
         GenerateFoldMatrix(pData(x)[,y],k,replicates,type,minclassize,
                           balanced, control)})
                           



















