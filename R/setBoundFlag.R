setBoundFlag <- function(X, cutoff) {
    x <- X[,3]
    x[x <= cutoff] <- 0
    x[x > cutoff] <- 1
    return(x)
}
