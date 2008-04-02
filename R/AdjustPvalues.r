### filename: AdjustPvalues.r
### Title: p-value adjustment by standard methods such as
###        those provided by the package 'multtest' or
###        the concept of 'qvalues'.
###
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 7.9.2007
### date(s) of updates:
#
### Brief description:
#
#  Takes a vector of pvalues and makes and adjustment
#  according to one of various methods.
#  -> mt.rawp2adjp(rawp, proc=c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY"))
#  qvalue.cal
#
#  Purpose: Finding a reasonale threshold.
#
#

#
### Further comments and notes:
#
#
#
###**************************************************************************###

AdjustPvalues <- function(pval,
                          method=c("BH", "qvalue", "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BY"))
{
 if(any(pval < 0 | pval > 1))
 stop("Raw pvalues are not from the unit interval \n")
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






 
 
 
