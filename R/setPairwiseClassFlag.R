setPairwiseClassFlag <- function(x, cutoff) {
    ws <- 0
    ws <- ifelse(x[,3] == 1 & x[,4] == 0 & (x[,1]-x[,2] >= cutoff),1,ws)
    ws <- ifelse(x[,3] == 0 & x[,4] == 1 & (x[,2]-x[,1] >= cutoff),2,ws)
    ws <- ifelse(x[,3] == 1 & x[,4] == 1,3,ws)
}
