probeAnnoFromESet <- function(eSet, probeLength) {

    stopifnot(inherits(eSet, "ExpressionSet"), is.numeric(probeLength))

    pos <- data.frame(fData(eSet)$CHR, fData(eSet)$PROBE_ID, fData(eSet)$START, probeLength)
    colnames(pos) <- c("CHROMOSOME", "PROBE_ID", "POSITION", "LENGTH")
    probeAnno <- posToProbeAnno(pos)
    return(probeAnno)
}
