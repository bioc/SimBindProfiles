setFlagIncreasedBinding <- function(x, cutoff) {
    ws <- 0
    ws <- ifelse(x[,3] == 1 & x[,4] == 1 & (x[,1]-x[,2] >= cutoff),1,ws)
}
