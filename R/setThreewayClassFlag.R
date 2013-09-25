setThreewayClassFlag <- function(x, cutoff) {
    ws <- 0
    ws <- ifelse(x[,4] == 1 & x[,5] == 0 & x[,6] == 0 & (x[,1]-x[,2] >= cutoff & x[,1]-x[,3] >= cutoff),1,ws)
    ws <- ifelse(x[,4] == 0 & x[,5] == 1 & x[,6] == 0 & (x[,2]-x[,1] >= cutoff & x[,2]-x[,3] >= cutoff),2,ws)
    ws <- ifelse(x[,4] == 0 & x[,5] == 0 & x[,6] == 1 & (x[,3]-x[,1] >= cutoff & x[,3]-x[,2] >= cutoff),3,ws)

    ws <- ifelse(x[,4] == 1 & x[,5] == 1 & x[,6] == 0 & (x[,1]-x[,3] >= cutoff & x[,2]-x[,3] >= cutoff),4,ws)
    ws <- ifelse(x[,4] == 0 & x[,5] == 1 & x[,6] == 1 & (x[,2]-x[,1] >= cutoff & x[,3]-x[,1] >= cutoff),5,ws)
    ws <- ifelse(x[,4] == 1 & x[,5] == 0 & x[,6] == 1 & (x[,1]-x[,2] >= cutoff & x[,3]-x[,2] >= cutoff),6,ws)

    ws <- ifelse(x[,4] == 1 & x[,5] == 1 & x[,6] == 1 ,7,ws)
}
