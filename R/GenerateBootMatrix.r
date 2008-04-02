### filename: GenerateBootMatrix.r
### Title: Genrate observation indices matrices for bootstrapping.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.8.2007
### date(s) of updates: 17.8.2007, 20.8.2007.
#
### Brief description:
#
#   Generates a matrix that is used for bootstrapping and
#   Repeated Ranking of Genes. Bootstrapping can be done
#   with various constraints (unique replications, minclassize, maxties etc.)
#   Up to now, algorithm is not determistic.
#   - Update from 20.8.2007 implements new balanced feature
#     (now two arguments: balancedclass, balancedsample)
#     balancedclass: Each class occurs the same time in each bootstrap
#                    replication.   
#     balancedsample: Each sample (array) occurs equally often (over all
#     replications).   
#
### Further comments and notes:
#
#   Very difficult and unperfect function.
#   Checks needed.
#   Can take a lot of time for big samplesize/k.
#   See also GenerateFoldMatrix.r.
#
###**************************************************************************###

### generic function

setGeneric("GenerateBootMatrix", function(x, y, replicates=50, type=c("unpaired", "paired", "onesample"),
            maxties=NULL, minclassize=2, balancedclass=FALSE, balancedsample=FALSE, control)
standardGeneric("GenerateBootMatrix"))

### signature: x=matrix, y=numeric.

setMethod("GenerateBootMatrix", signature(x="matrix", y="numeric"),
function(x, y, replicates=50, type=c("unpaired", "paired", "onesample"),
         maxties=NA, minclassize=2, balancedclass=FALSE, balancedsample=FALSE, 
         control){
  ly <- length(y)
  if(missing(control)) control <- samplingcontrol(candreplicates=replicates*3)
  candreplicates <- with(control, candreplicates)
  maxiter <- with(control, maxiter)
  UNIQUE <- FALSE
  iter <- 0
  type <- match.arg(type)
  if( !is.element(type, eval(formals(GenerateBootMatrix)$type)))
  stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'.")
  taby <- table(y)
  if(length(taby) != 2 & type != "onesample")
  stop("y has not exactly two levels, but 'type' is not 'onesample'.")
  if(length(taby) != 1 & type == "onesample")
  warning("y has not exactly one level, but 'type' is 'onesample'. ")
  if(!balancedsample){
      if(type=="paired") samp <- ly/2
      else samp <- ly
      REPLACE <- TRUE
  }
  else{ 
    if(type == "paired") samp <- rep(1:(ly/2), times=replicates)
    else samp <- rep(1:ly, times=replicates)
     if(balancedclass | minclassize >= 2 | !is.na(maxties)){
      warning("Further restrictions currently not admitted if
               balancedsample = 'TRUE'. \n")
        balancedclass <- FALSE
        minclassize <- 1       
        maxties <- NA       
      }
      candreplicates <- 1
      REPLACE <- FALSE
    }
  if(type == "paired"){
  if(length(unique(y[1:taby[1]])) != 1 | length(unique(y[-c(1:taby[1])])) != 1)
           stop("Incorrect coding for type='paired'. \n")
  }
  if(!is.na(maxties)){
  maxties <- as.integer(maxties)
  if(maxties < 1) stop("Invalid value specified for argument 'maxties'.")
  }
  if(type == "unpaired"){
    if(balancedclass) minclassize <- floor(ly/2)
        if(minclassize < 1){
         minclassize <- 1
         warning("Invalid minclassize for type='unpaired'. Reset to 1")
        }
    if(is.na(maxties) || maxties >= (ly-minclassize-1)){
        checkfun <- function(z){
         row <- y[z]
         tab <- table(row)
         if(tab[1] >= minclassize && tab[2] >= minclassize) return(z)
         else return(rep(NA,ly))}
       }
    else{
     checkfun <- function(z){
         row <- y[z]
         tab <- table(row)
         tab2 <- table(z)
         if(tab[1] >= minclassize && tab[2] >= minclassize && (max(tab2) <= (maxties+1))) return(z)
         else return(rep(NA,ly))}
     }
      while(iter < maxiter && !UNIQUE){
      bootm <- replicate(candreplicates, sample(samp, replace=REPLACE))
      if(balancedsample) bootm <- matrix(bootm, nrow=ly, ncol=replicates)
      iter <- iter+1
      bootm <- apply(bootm, 2, checkfun)
      indbootm <- apply(bootm, 2, function(z) any(!is.na(z)))
      if(length(indbootm) < replicates){
        if(iter < maxiter){
        cat("Not enough appropriate replications foud. Repeating... \n")
        next
        }
        else{
        warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
        uniqbootm <- unique(bootm[,indbootm], MARGIN=2)
        replicates <- ncol(uniqbootm)
        }
        }
      else{
      bootm <- bootm[,indbootm]
      uniqbootm <- unique(bootm, MARGIN=2)
      UNIQUE <- (ncol(uniqbootm) >= replicates)
      if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
          replicates <- ncol(bootm)
          }
         }
        }
       }
       bootm <- uniqbootm
       if(ncol(bootm)>replicates) bootm <- bootm[,sample(ncol(bootm), replicates)]
      }
  if(type == "onesample"){
   if(is.na(maxties) || maxties >= (ly-1)){
       while(iter < maxiter && !UNIQUE){
        bootm <- replicate(candreplicates, sample(samp, replace=REPLACE))
        if(balancedsample) bootm <- matrix(bootm, nrow=ly, ncol=replicates) 
        iter <- iter+1
        uniqbootm <- unique(bootm, MARGIN=2)
        UNIQUE <- (ncol(uniqbootm) >= replicates)
        if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
          bootm <- uniqbootm
          replicates <- ncol(uniqbootm)
          }
          }
         else bootm <- uniqbootm[,sample(ncol(uniqbootm), replicates)]
        }
       }
   else{
    checkfun <- function(z){ row <- rep(NA, ly); tab <-  table(z)
                              if(max(tab) > (maxties+1)) return(row)
                              else return(z)}
      while(iter < maxiter && !UNIQUE){
      bootm <- replicate(candreplicates, sample(samp, replace=REPLACE))
      if(balancedsample) bootm <- matrix(bootm, nrow=ly, ncol=replicates)
      iter <- iter+1
      bootm <- apply(bootm, 2, checkfun)
      indbootm <- apply(bootm, 2, function(z) any(!is.na(z)))
      if(length(indbootm) < replicates){
        if(iter < maxiter){
        cat("Not enough appropriate replications foud. Repeating... \n")
        next
        }
        else{
        warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
        uniqbootm <- unique(bootm[,indbootm], MARGIN=2)
        replicates <- ncol(uniqbootm)
        }
        }
      else{
      bootm <- bootm[,indbootm]
      uniqbootm <- unique(bootm, MARGIN=2)
      UNIQUE <- (ncol(uniqbootm) >= replicates)
      if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
          replicates <- ncol(bootm)
          }
         }
        }
       }
       bootm <- uniqbootm
       if(ncol(bootm)>replicates) bootm <- bootm[,sample(ncol(bootm), replicates)]
      }
      balancedclass <- FALSE
      minclassize <- 0
     }
  if(type == "paired"){
   if(is.na(maxties) || maxties >= (ly/2-1)){
       while(iter < maxiter && !UNIQUE){
        bootm <- replicate(candreplicates, sample(samp, replace=REPLACE))
        if(balancedsample) bootm <- matrix(bootm, nrow=ly/2, ncol=replicates)
        iter <- iter+1
        uniqbootm <- unique(bootm, MARGIN=2)
        UNIQUE <- (ncol(uniqbootm) >= replicates)
        if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
          bootm <- uniqbootm
          replicates <- ncol(uniqbootm)
          }
          }
         else bootm <- uniqbootm[,sample(ncol(uniqbootm), replicates)]
        }
       }
   else{
    checkfun <- function(z){ row <- rep(NA, ly/2); tab <-  table(z)
                              if(max(tab) > (maxties+1)) return(row)
                              else return(z)}
      while(iter < maxiter && !UNIQUE){
      bootm <- replicate(candreplicates, sample(samp, replace=REPLACE))
      if(balancedsample) bootm <- matrix(bootm, nrow=ly/2, ncol=replicates)
      iter <- iter+1
      bootm <- apply(bootm, 2, checkfun)
      indbootm <- apply(bootm, 2, function(z) any(!is.na(z)))
      if(length(indbootm) < replicates){
        if(iter < maxiter){
        cat("Not enough appropriate replications foud. Repeating... \n")
        next
        }
        else{
        warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
        uniqbootm <- unique(bootm[,indbootm], MARGIN=2)
        replicates <- ncol(uniqbootm)
        }
        }
      else{
      bootm <- bootm[,indbootm]
      uniqbootm <- unique(bootm, MARGIN=2)
      UNIQUE <- (ncol(bootm) >= replicates)
      if(!UNIQUE){
          if(iter < maxiter) cat("Not enough appropriate replications foud. Repeating... \n")
          else{
          warning("Desired Number of replications could not be generated. \n Check all arguments for validity.")
          replicates <- ncol(bootm)
          }
         }
        }
       }
       bootm <- uniqbootm
       if(ncol(bootm)>replicates) bootm <- bootm[,sample(ncol(bootm), replicates)]
      }
      bootm <- rbind(bootm, bootm+ly/2)
      balancedclass <- TRUE
      minclassize <- ly/2
     }
    new("BootMatrix", bootmatrix=bootm, replicates=as.integer(replicates), type=type,
       maxties=maxties, minclassize=as.integer(minclassize), balancedclass=balancedclass,
        balancedsample = balancedsample)
}
)

### signature: x=matrix, y=factor.

setMethod("GenerateBootMatrix", signature(x="matrix", y="factor"),
function(x, y, replicates=50, type=c("unpaired", "paired", "onesample"),
         maxties=NA, minclassize=2, balancedclass=FALSE, balancedsample=FALSE, control){
         GenerateBootMatrix(x,as.numeric(y), replicates, type, maxties, minclassize, 
                            balancedsample, balancedclass, control)})

### signature: x=ExpressionSet, y=character

setMethod("GenerateBootMatrix", signature(x="ExpressionSet", y="character"),
function(x, y, replicates=50, type=c("unpaired", "paired", "onesample"),
         maxties=NA, minclassize=2, balancedsample=FALSE, balancedclass=FALSE, 
         control){
         GenerateBootMatrix(exprs(x),pData(x)[,y],replicates,type,maxties,minclassize,
                           balancedclass,balancedsample,control)})
                           






