readSgrFiles <- function(X, dataPath=getwd(), fileExt=".txt", normalise=TRUE) {

    signal <- c()
    rowCheck <- NULL
    chromosomeCheck <- NULL
    for(i in 1:length(X)) {
        dataFile <- paste(dataPath, "/", X[i], fileExt, sep="")
        stopifnot(file.exists(dataFile))
        cat("Reading", dataFile, "...\n")
        sgr <- read.delim(dataFile , as.is=TRUE, header=FALSE)
        order.sgr <- sgr[order(sgr[,1], sgr[,2]),]
        colnames(order.sgr) <- c("CHR","START",X[i])
       
        if((!is.null(rowCheck)) && (rowCheck != nrow(sgr)))
          stop("Files do not contain the same number of rows!\n")
        if((!is.null(chromosomeCheck)) && length(setdiff(sort(unique(order.sgr[,1])), chromosomeCheck))  > 0) {
          stop("Chromosome names don't match between files!\n",
               "Chromsome names in first file:\n", paste(chromosomeCheck, "", sep=" "), "\n",
               "Found in", dataFile ,":\n", paste(setdiff(sort(unique(order.sgr[,1])), chromosomeCheck), "", sep=" "), "\n")
        }
        if(is.null(rowCheck)) rowCheck <- nrow(signal)
        if(is.null(chromosomeCheck)) chromosomeCheck <- sort(unique(order.sgr[,1]))
        signal <- cbind(signal, order.sgr[,3])
    }
  
    PROBE_ID <- paste(order.sgr[,1], order.sgr[,2], sep="")
    probe.location <- data.frame(PROBE_ID, order.sgr[,1], order.sgr[,2], row.names = PROBE_ID)
    colnames(probe.location) <- c("PROBE_ID", "CHR", "START")
    targets <- data.frame(X, "dataType", row.names=X)
    colnames(targets) <- c("FileName", "tiling array")
    if(normalise == TRUE) {
        cat("Performing quantile normalisation...\n")
        signal <- normalizeBetweenArrays(as.matrix(signal), method="quantile")
    }
    cat("Building ExpressionSet...\n\n")
    myPD <- new("AnnotatedDataFrame", data = targets, varMetadata = data.frame(varLabel = colnames(targets),
        row.names = colnames(targets)))
    myEset <- new("ExpressionSet", exprs = signal, phenoData = myPD)
    featureNames(myEset) <- make.names(probe.location[["PROBE_ID"]], unique = TRUE)
    featureData(myEset) <- as(probe.location, "AnnotatedDataFrame")
	if (validObject(myEset)) {
		return(myEset)
	}
}
