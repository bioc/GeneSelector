### filename: AggregateSimple.r
### Title: Aggregation using simple statistics.
### Author: Martin Slawski
### email: <Martin.Slawski@campus.lmu.de>
### date of creation: 31.8.2007
### date(s) of updates: 3.9.2007
### major revision:  3.12.2008
#
### Brief description:
#
#  Aggregation is based on simple (weighted) descriptive
#  statistics; weights depend on stability measures.
#
#
### Further comments and notes:
#
#
###**************************************************************************###

setGeneric("AggregateSimple", function(RR, measure = c("mode", "mean", "trimmed.mean", "median", "quantile"), q=NULL, trim = NULL)
            standardGeneric("AggregateSimple"))
            

setMethod("AggregateSimple", signature(RR="RepeatedRanking"),
            function(RR, measure = c("mode", "mean", "trimmed.mean", "median", "quantile"), q=NULL, trim = NULL){

  R <- cbind(RR@original@ranking, RR@rankings) ### note: completely symmetric treatment.
  p <- nrow(R)

  measure <- match.arg(measure)
  if(!is.element(measure, eval(formals(AggregateSimple)$measure)))
  stop("Invalid 'measure' specified \n")

  if(is.null(q) & measure  == "quantile"){
    measure <- "median"
    warning("'measure' is 'quantile', but argument 'q' is NULL; set to 0.5 \n")
  }
  
  if(!is.null(trim) && (trim < 0 || trim > 0.5))
   stop("'trim' may range from 0 to 0.5 \n")
  
  if(is.null(trim) & measure  == "trimmed.mean"){
    measure <- "mean"
    warning("'measure' is 'trimmed.mean', but argument 'trim' is NULL; set measure to 'mean' \n")
  }

  aggregatefun <-  switch(measure,  mode=function(rb){
                                        tab <- table(rb)
                                        tabnames <- names(tab)
                                        maxtab <- tab[which.max(tab)]
                                        modes <- as.numeric(tabnames[tab == maxtab])
                                        lmodes <- length(modes)
                                        if(lmodes == 1) return(modes)
                                        else return(modes[which.min(abs(modes - rb[1]))])
                                        },
                                        mean=function(rb) mean(rb),
                                        trimmed.mean = function(rb) mean(rb, trim = trim),
                                        median=function(rb) median(rb),
                                        quantile=function(rb) quantile(rb, q))


  summary <- apply(R, 1, aggregatefun)
  

  measureout <- switch(measure, mode = "most frequent",
                                mean = "mean",
                                trimmed.mean = paste(trim, "trimmed mean", sep="%-"),
                                median = "median",
                                quantile = paste(q, "quantile", sep = "-"))
  ### note: ties 'ignored'.
  new("AggregatedRanking", ranking = rank(summary, ties.method = "first"), type= "simple", measure = measureout, method = RR@original@method)
 }
)

