### filename: GetAlpha.r
### Title:
### Determine weight decay parameter for stability scores
### (Lottaz, Linear Model with weights)
###
###
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 10.9.2007
### date(s) of updates:
#
### Brief description:
#
#   Finds an appropriate value for weight parameter alpha
#   based on a nonlinear regression of p-values on ranks.
#
#

#
### Further comments and notes:
#
#   Pval is assumed to correspond to ranking.
#   Input can either be from a ranking object (indices, sorted pvalues)
#   or a list of ranks with correspong p-values or adjusted p-values,
#   respectively.
#
###**************************************************************************###

GetAlpha <- function(ranking, pval, alpha0=0.01){
if(length(ranking) != length(pval))
stop("length of 'ranking and 'pval' differ \n")
p <- length(ranking)
if(any(pval < 0 | pval > 1))
stop("pvalues must between zero and one")
x <- match(1:p, ranking)
y <- pval[x]
nlsm <- nls(y ~ 1-exp(-alpha*x), start=list(alpha=0.1))
return(coef(nlsm))
}

