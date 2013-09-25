increasedBindingRegions <- function(xSet, sgrset=c(1,2), bound.cutoff, diff.cutoff, 
                            probes, probe.max.spacing, writeBedFile=TRUE) {

    stopifnot(inherits(xSet, "ExpressionSet"), is.numeric(bound.cutoff), is.numeric(diff.cutoff),
        is.numeric(probes), is.numeric(probe.max.spacing))

    if(length(sgrset) != 2) stop("You must specify 2 data sets from the xSet object!")
    for(i in 1:length(sgrset)) {
        if(is.na(sampleNames(xSet)[sgrset[i]])) 
          stop("Your selected sgrset ",sgrset[i]," doesn't exist in the xSet object!\n")
    }

    X <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[1]])
    Y <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[2]])
    name1 <- sampleNames(xSet)[sgrset[1]]
    name2 <- sampleNames(xSet)[sgrset[2]]
    cat("\nIncreased Binding comparison of",name1, "vs", name2, "\n\n")
    table <- data.frame(X[,3], Y[,3], setBoundFlag(X, cutoff=bound.cutoff), 
                        setBoundFlag(Y, cutoff=bound.cutoff))
    class.table <- data.frame(as.character(X[,1]), X[,2], setFlagIncreasedBinding(table, cutoff=diff.cutoff), 
                              X[,3], Y[,3], stringsAsFactors = FALSE)
    colnames(class.table) <- c("chr", "start", "class", name1, name2)

    if(nrow(class.table[class.table$class == 1,]) == 0) {
        cat("Found no probes in class 1\n")
    }
    if(nrow(class.table[class.table$class == 1,]) > 1) {
        cat("Filter data into regions...\n")
        data.class <- filterRegionsPairwise(class.table, class.group=1, probes=probes, 
                                            probe.max.spacing=probe.max.spacing)
        bedname <- paste(name1,".vs.",name2, ".increasedBinding_b", bound.cutoff, 
                         "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
        if(!is.null(data.class)) {
            if(writeBedFile) {
                bed <- cbind(data.class[,1:3], paste(name1, ".increasedBinding", sep=""), 
                             data.class[,5])
                bed <- bed[order(bed[,5], decreasing = TRUE), ]
                cat("Writing",bedname, ",regions =",nrow(bed),"...\n")
                write.table(bed ,file=bedname, col.names=F, sep="\t", row.names=F, quote=F)
            }
            boundRegion <- cbind(paste(name1, ".increasedBinding", sep=""), 1, 
                                 data.class[,1:3], data.class[,5:6])
            colnames(boundRegion) <- c("name", "class.group", "chr", "start", "end", "score", "nProbes")
        }
        return(boundRegion)
    }
}
