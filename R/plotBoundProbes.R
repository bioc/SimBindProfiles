plotBoundProbes <- function(xSet, sgrset=c(1,2), method=c("pairwise","compensation","increasedBinding"),
        bound.cutoff, diff.cutoff, cols=NULL, pcex=2) {

    stopifnot(inherits(xSet, "ExpressionSet"), is.numeric(bound.cutoff), is.numeric(diff.cutoff))

    if(length(sgrset) < 2) stop("You must specify at least 2 data sets from the xSet object!")   
    for(i in 1:length(sgrset)) {
        if(is.na(sampleNames(xSet)[sgrset[i]])) 
          stop("Your selected sgrset ",sgrset[i]," doesn't exist in the xSet object!\n")
    }

    method = match.arg (method)
    X <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[1]])
    Y <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[2]])
    Z <- cbind(fData(xSet)[2:3], exprs(xSet)[,sgrset[3]])
    name1 <- sampleNames(xSet)[sgrset[1]]
    name2 <- sampleNames(xSet)[sgrset[2]]
    name3 <- sampleNames(xSet)[sgrset[3]]

    if(method == "pairwise") {
        bound <- cbind(setBoundFlag(X, cutoff=bound.cutoff), setBoundFlag(Y, cutoff=bound.cutoff))
        table <- cbind(X[,3], Y[,3], bound)
        class <- setPairwiseClassFlag(table, cutoff=diff.cutoff)
        if(is.null(cols) | (length(cols) != 3)) {
            cols <- c("red", "green", "grey")
        }
        plot(x=X[,3], y=Y[,3], pch=".", xlab=name1, ylab=name2, main="DIFF BOUND")
        points(x=X[class==1,3], y=Y[class==1,3], pch=".", cex=pcex, col=cols[1])
        points(x=X[class==2,3], y=Y[class==2,3], pch=".", cex=pcex, col=cols[2])
        points(x=X[class==3,3], y=Y[class==3,3], pch=".", cex=pcex, col=cols[3])
    }

    if(method == "increasedBinding") {
        bound <- cbind(setBoundFlag(X, cutoff=bound.cutoff), setBoundFlag(Y, cutoff=bound.cutoff))
        table <- cbind(X[,3], Y[,3], bound)
        class <- setFlagIncreasedBinding(table, cutoff=diff.cutoff)
        if(is.null(cols)) {
            cols <- "blue"
        }
        plot(x=X[,3], y=Y[,3], pch=".", xlab=name1, ylab=name2, main="INCREASED BINDING")
        points(x=X[class==1,3], y=Y[class==1,3], pch=".", cex=pcex, col=cols[1])
    }

    if(method == "compensation") {
        if(length(sgrset) != 3) stop("You must select 3 data sets for compensation method!\n")
        bound <- cbind(setBoundFlag(X, cutoff=bound.cutoff), setBoundFlag(Y, cutoff=bound.cutoff), setBoundFlag(Z, cutoff=bound.cutoff))
        table <- cbind(X[,3], Y[,3], Z[,3], bound)
        class <- setFlagCompensation(table, cutoff=diff.cutoff)
        avgXY <- (X[,3] + Y[,3])/2
        if(is.null(cols)) {
            cols <- "orange"
        }
        plot(x=avgXY, y=Z[,3], pch=".", xlab=paste("average ", name1, ".", name2, sep=""), ylab=name3, main="COMPENSATION")
        points(x=avgXY[class==1], y=Z[class==1,3], pch=".", cex=pcex, col=cols[1])
    }
}
