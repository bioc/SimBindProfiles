setFlagCompensation <- function(x, cutoff) {
    ws <- 0
    ws <- ifelse(x[,4] == 1 & x[,5] == 1 & x[,6] == 0 & (((x[,1]+x[,2])/2 - x[,3]) >= cutoff),1,ws) # 1,1,0
}
