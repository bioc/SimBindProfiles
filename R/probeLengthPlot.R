probeLengthPlot <- function(xSet, sgrset=1, chr=NULL, bound.cutoff, probe.max.spacing=200, xlim.max=25) {

    stopifnot(inherits(xSet, "ExpressionSet"), is.numeric(bound.cutoff))
    if(is.na(sampleNames(xSet)[sgrset])) stop("Your selected sgrset doesn't exist in the xSet object!\n")

    X <- data.frame(as.character(fData(xSet)[,2]), fData(xSet)[,3], exprs(xSet)[,sgrset],
                    stringsAsFactors = FALSE)
    colnames(X) <- c("CHR", "START", "Ratio")
    if(!is.null(chr)) {
      chr <- sub("chr", "", chr)
      chr <- paste("chr", chr, sep="")
      match.arg(chr, unique(fData(xSet)$CHR))
      X <- subset(X, X$CHR == chr)
    }
    c <- X[,3]
    
    pick.cutoff = which(c >= bound.cutoff)
    peaks.raw = list()
    p = c()
    for(i in 1:length(pick.cutoff)){
        inx = pick.cutoff[i]
        if(length(p)==0 || (tail(p,1) == inx-1 && (X[inx,2]-X[inx-1,2]) <= probe.max.spacing)){
            p = c(p, inx)
        } else {
            if(length(p) !=0){
                peaks.raw[[length(peaks.raw)+1]] = p
                p = c(inx)
            }
        }
    }
    if(length(p) !=0) {
        peaks.raw[[length(peaks.raw)+1]] = p
    }
    nProbes <- c()
    for(i in 1:length(peaks.raw)){
        p.nProbes = length(peaks.raw[[i]])
        nProbes[i] = p.nProbes
    }
    ptable<- table(nProbes)
    plot(ptable, xlim=c(1,xlim.max), main=paste(sampleNames(xSet)[sgrset], "with bound.cutoff =", 
        bound.cutoff, sep=" "), lwd=5, ylab="Frequency")
}
