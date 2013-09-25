eSetScatterPlot <- function(xSet) {
    stopifnot(inherits(xSet, "ExpressionSet"))
    X <- exprs(xSet)

    upPan <- function(...){
        par(new=TRUE);smoothScatter(..., nrpoints=0,)
    }
    lowPan <- function(x, y, ...) {
        text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
        signif(cor(x, y), 2), cex = 1.5)
    }
    pairs(X[, 1:(ncol(X))], pch = ".", lower.panel = lowPan, 
        upper.panel = upPan)
}
