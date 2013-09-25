threewayRegions <- function(xSet, sgrset=c(1,2,3), bound.cutoff, diff.cutoff, probes, 
                            probe.max.spacing, writeBedFile=TRUE) {

    stopifnot(inherits(xSet, "ExpressionSet"), is.numeric(bound.cutoff), is.numeric(diff.cutoff),
        is.numeric(probes), is.numeric(probe.max.spacing))

    if(length(sgrset) != 3) stop("You must specify 3 data sets from the xSet object!")
    for(i in 1:length(sgrset)) {
        if(is.na(sampleNames(xSet)[sgrset[i]])) 
          stop("Your selected sgrset ",sgrset[i]," doesn't exist in the xSet object!\n")
    }

    X <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[1]])
    Y <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[2]])
    Z <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[3]])
    name1 <- sampleNames(xSet)[sgrset[1]]
    name2 <- sampleNames(xSet)[sgrset[2]]
    name3 <- sampleNames(xSet)[sgrset[3]]
    cat("\nThree way comparison of",name1, "vs", name2, "vs", name3, "\n\n")
    table <- data.frame(X[,3], Y[,3], Z[,3], setBoundFlag(X, cutoff=bound.cutoff), 
                        setBoundFlag(Y, cutoff=bound.cutoff), setBoundFlag(Z, cutoff=bound.cutoff))
    class.table <- data.frame(as.character(X[,1]), X[,2], setThreewayClassFlag(table, cutoff=diff.cutoff), 
                              X[,3], Y[,3], Z[,3], stringsAsFactors = FALSE)
    colnames(class.table) <- c("chr", "start", "class", name1, name2, name3)
    cat("Filter data into regions...\n")
    boundRegionSummary <- c()

    for (i in 1:7) {
        if(nrow(class.table[class.table$class == i,]) == 0) {
            cat("Found no probes in class ", i,"\n")
        }
        if(nrow(class.table[class.table$class == i,]) > 1) {
            data.class <- filterRegionsThreeway(class.table, class.group=i, 
                                                probes=probes, probe.max.spacing=probe.max.spacing)
            if(i == 1) {
                bedname <- paste(name1,".vs.",name2, ".", name3, ".unique_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name1, ".unique", sep="")
            }
            if(i == 2) {
                bedname <- paste(name2,".vs.",name1, ".", name3, ".unique_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name2, ".unique", sep="")
            }
            if(i == 3) {
                bedname <- paste(name3, ".vs.", name1, ".", name2, ".unique_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name3, ".unique", sep="")
            }
            if(i == 4) {
                bedname <- paste(name1,".",name2, ".vs.", name3, ".common_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name1, ".", name2, ".common", sep="")
            }
            if(i == 5) {
                bedname <- paste(name2,".",name3, ".vs.", name1, ".common_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name2, ".", name3, ".common", sep="")
            }
            if(i == 6) {
                bedname <- paste(name1,".",name3, ".vs.", name2, ".common_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name1, ".", name3, ".common", sep="")
            }
            if(i == 7) {
                bedname <- paste(name1,".",name2, ".", name3, ".common_b", bound.cutoff, 
                                 "d",diff.cutoff, "v",probes,"g",probe.max.spacing,".bed", sep="")
                bname <- paste(name1, ".", name2, ".", name3, ".common", sep="")
            }
            if(!is.null(data.class)) {
                if(writeBedFile) {
                    bed <- cbind(data.class[,1:3], bname, data.class[,5])
                    bed <- bed[order(bed[,5], decreasing = TRUE), ]
                    cat("Writing",bedname, ",regions =",nrow(bed),"...\n")
                    write.table(bed ,file=bedname, col.names=F, sep="\t", row.names=F, quote=F)
                }
                boundRegion <- cbind(bname, i, data.class[,1:3], data.class[,5:6])
                boundRegionSummary <- rbind(boundRegionSummary, boundRegion) 
            }
        }
    }
    colnames(boundRegionSummary) <- c("name", "class.group", "chr", "start", "end", "score", "nProbes")
    return(boundRegionSummary)
}
