### filename: Internals.r
### Title: Internal functions not to be called directly by the user.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 16.8.2007
### date(s) of updates: 27.8.2007; 31.8.2007; 3.9.2007; 6.9.2007; 10.9.2007
#
### Brief description:
#
#   Helper functions for the important methods.
#
### Further comments and notes:
#
#
#
###**************************************************************************###

### [1]

combn <- function (x, m)
{
    if (length(m) > 1) {
        warning(paste("Argument m has", length(m), "elements: only the first used"))
        m <- m[1]
    }
    if (m < 0)
        stop("m < 0")
    if (m == 0)
        return(vector(mode(x), 0))
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) ==
        x)
        x <- seq(x)
    n <- length(x)
    if (n < m)
        stop("n < m")
    e <- 0
    h <- m
    a <- 1:m
    count <- choose(n,m)
    out <- vector("list", count)
    out[[1]] <- x[a]
    dim.use <- NULL
    if (count > 1) dim.use <- c(m, count)
    i <- 2
    nmmp1 <- n - m + 1
    mp1 <- m + 1
    while (a[1] != nmmp1) {
        if (e < n - h) {
            h <- 1
            e <- a[m]
            j <- 1
        }
        else {
            h <- h + 1
            e <- a[mp1 - h]
            j <- 1:h
        }
        a[m - h + j] <- e + j
        out[[i]] <- x[a]
        i <- i + 1
    }
        if (is.null(dim.use))
            out <- unlist(out)
        else out <- array(unlist(out), dim.use)
    return(out)
}

### [2]

samplingcontrol <- function(candreplicates, maxiter=5){
list(candreplicates=candreplicates, maxiter=maxiter)}

### [3]  function for the Overlap Score of Lottaz

overlap <- function (x1, x2, n) 
{
    r <- match(x1[1:n], x2[1:n])
    overlapRanks <- pmax(r, 1:n)
    tmp <- table(overlapRanks)
    x <- integer(n)
    x[as.integer(names(tmp))] <- tmp
    x <- cumsum(x)
    return(x)
}

### [4] function for union count

unioncount <- function(M){
     p <- nrow(M)
     B <- ncol(M)
     active <- unique(M[1,])
     count <-  numeric(p)
     count[1] <- length(active)
      for(j in 2:p){
       mj <- M[j,] %in% active
       if(length(mj) == 0) active <- c(active, M[j,])
       else active <- c(active, unique(M[j,!mj]))
       count[j] <- length(active)
       if(count[j]  >=  p){
        count[j] <- p
        break
        }

     }
     if((j+1) <= p)  count[(j+1):p] <- p
       count
    }
    
### [5] function for p-value adjustment (for GeneSelector)

AdjustPvalues <- function(pval,
                          method=c("BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"))
{
 if(any(pval < 0 | pval > 1))
 stop("Raw p-values are not contained in the unit interval \n")
 method <- match.arg(method)
 if(!is.element(method, eval(formals(AdjustPvalues)$method)))
 stop("Invalid method specified. \n")
 if(method != "qvalue"){
 require(multtest, quietly=TRUE)
 outp <- mt.rawp2adjp(pval, proc=method)
 adjpval <- outp$adjp[order(outp$index),-1]
 }
 else{
   require(siggenes, quietly=TRUE)
   p0 <- pi0.est(pval)
   adjpval <- siggenes:::qvalue.cal(pval, p0$p0, version = 2)
   }
 return(adjpval)
}



### [4] function for posterior distributions for AggregateBayes.r

#posterior1 <- function(x, omega, tau, p){
#lx <- length(x)
#rb <- numeric(lx-1)
#r0 <- x[1]
#rb <- x[2:lx]
#mu0 <- weighted.mean(rb, omega)
#omega <- omega/mean(omega)
#sigmahat <- mad(abs(rb-mu0)*omega, constant=1)
#
#if(sigmahat < .Machine$double.eps){
#    post <- 1
#    names(post) <- mu0
#    return(post)
#   }
#else{
#lower <- max(1, r0-round(3*tau))
#upper <- min(r0+round(3*tau), p)
#supp <- lower:upper
#prior <- dnorm(supp, mean=r0, sd=tau)
#lik <- exp(sapply(supp, function(zz)
#                  sum(dnorm(rb[rb %in% supp], mean=zz, sd=sigmahat, log=TRUE))))
#unpost <- lik*prior
#if(all(unpost < .Machine$double.eps)) unpost <- prior
#post <- unpost/sum(unpost)
#names(post) <- as.character(supp)
#return(post)
#}
#}


#posterior2 <- function(x, omega, p){
#lx <- length(x)
#tau <- x[1]
#rb <- numeric(lx-2)
#r0 <- x[2]
#rb <- x[3:lx]
#mu0 <- weighted.mean(rb, omega)
#omega <- omega/mean(omega)
#sigmahat <- mad(abs(rb-mu0)*omega, constant=1)

#if(sigmahat < .Machine$double.eps){
#    post <- 1
#    names(post) <- mu0
#    return(post)
#   }
#else{
#lower <- max(1, r0-round(3*tau))
#upper <- min(r0+round(3*tau), p)
#supp <- lower:upper
#prior <- dnorm(supp, mean=r0, sd=tau)
#lik <- exp(sapply(supp, function(zz)
#                  sum(dnorm(rb[rb %in% supp], mean=zz, sd=sigmahat, log=TRUE))))
#unpost <- lik*prior
#if(all(unpost < .Machine$double.eps)) unpost <- prior
#post <- unpost/sum(unpost)
#names(post) <- as.character(supp)
#return(post)
#}
#}

#posterior3 <- function(x, omega, tau, p){
#lx <- length(x)
#rb <- numeric(lx-1)
#r0 <- x[1]
#rb <- x[2:lx]
#mu0 <- weighted.mean(rb, omega)
#omega <- omega/mean(omega)
#sigmahat <- sqrt(sum(omega*(rb - mu0)^2)/sum(omega))

#if(sigmahat < .Machine$double.eps){
#    post <- 1
#    names(post) <- mu0
#    return(post)
#   }
#else{
#lower <- max(1, r0-round(3*tau))
#upper <- min(r0+round(3*tau), p)
#supp <- lower:upper
#prior <- dnorm(supp, mean=r0, sd=tau)
#lik <- exp(sapply(supp, function(zz)
#                  sum(dnorm(rb[rb %in% supp], mean=zz, sd=sigmahat, log=TRUE))))
#unpost <- lik*prior
#if(all(unpost < .Machine$double.eps)) unpost <- prior
#post <- unpost/sum(unpost)
#names(post) <- as.character(supp)
#return(post)
#}
#}

#posterior4 <- function(x, omega, p){
#lx <- length(x)
#tau <- x[1]
#rb <- numeric(lx-2)
#r0 <- x[2]
#rb <- x[3:lx]
#mu0 <- weighted.mean(rb, omega)
#omega <- omega/mean(omega)
#sigmahat <- sqrt(sum(omega*(rb - mu0)^2)/sum(omega))

#if(sigmahat < .Machine$double.eps){
#    post <- 1
#    names(post) <- mu0
#    return(post)
#   }
#else{
#lower <- max(1, r0-round(3*tau))
#upper <- min(r0+round(3*tau), p)
#supp <- lower:upper
#prior <- dnorm(supp, mean=r0, sd=tau)
#lik <- exp(sapply(supp, function(zz)
#                  sum(dnorm(rb[rb %in% supp], mean=zz, sd=sigmahat, log=TRUE))))
#unpost <- lik*prior
#if(all(unpost < .Machine$double.eps)) unpost <- prior
#post <- unpost/sum(unpost)
#names(post) <- as.character(supp)
#return(post)
#}
#}

### [5] Tuning for the SoftThresholded t-Statistic from Wu.
###     adapted from K. Strimmer. 

getLambda <- function(di, si, method = c("cor", "lowess"))
{
    if (method == "lowess") {
        Lambda <- c(0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50)
        rL <- sapply(Lambda, function(lambda) {
            mu12 <- ifelse(abs(di) > lambda, di - sign(di) * lambda, 
                0)
            tmp <- mu12/sqrt(2 * si^2 + lambda^2/3)
            a <- lowess(si, tmp, f = 2/3)$y
            mean((tmp[order(si)] - a)^2)/var(tmp)
        })
        Lam.Opt = Lambda[order(-rL)[1]]
    }
    if (method == "cor") {
        Lambda <- seq(0, 100, length = 100)
        crL <- sapply(Lambda, function(lambda) {
            mu12 <- ifelse(abs(di) > lambda, di - sign(di) * lambda, 
                0)
            tmp <- mu12/sqrt(2 * si^2 + lambda^2/3)
            idm <- order(-si)[1]
            return(cor(tmp[-idm], si[-idm]))
        })
        Lam.Opt <- Lambda[order(abs(crL))[1]]
        }
    
    return(Lam.Opt)
}

#### [6] For the GeneSelector Distance Plot

#characterplot <- function(char, x, y, deltax, deltay, cex=1){
#spltchar <- unlist(strsplit(char, ""))
#for(s in seq(along=spltchar)) points(x+(s-1)*deltax, y+deltay, pch=spltchar[s], cex=cex)
#}
