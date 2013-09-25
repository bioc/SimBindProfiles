findBoundCutoff <- function(xSet, method=c("normalNull","twoGaussiansNull"), mean.method="mode", 
                    pvalue=FALSE, fdr=FALSE, pvalPlot=FALSE) {

    stopifnot(inherits(xSet, "ExpressionSet"))
    signal <- exprs(xSet)[,1]
    if(method == "normalNull") {
        pval <- normalNull(signal, mean.method = mean.method)
    }
    if(method == "twoGaussiansNull") {
         pval <- twoGaussiansNull.BF(signal)
    }
    if(pvalPlot) {
        hist(pval, breaks=1000, freq=FALSE, border="grey", main="Histogram of p", 
            xlab="pvalue")
    }
    if(pvalue) {
        cut <- min(signal[pval <= pvalue], na.rm=TRUE)
    }
    if(fdr) {
        adj.pval <- p.adjust(pval, method = "fdr")
        cut <- min(signal[adj.pval <= fdr], na.rm=TRUE)
    }
    cutoff <- round(cut,2)
    cat("Using bound.cutoff =", cutoff, "\n")
    return(cutoff)
}
