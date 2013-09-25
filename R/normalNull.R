normalNull = function(x, mean.method=c("zero", "mode")) {

    x <- as.numeric(na.omit(x))
    mean.method = match.arg (mean.method)
    if (mean.method == "zero") {
        x.mean = 0
        }
    if (mean.method == "mode") {
        dd <- density(x)
        x.mean = dd$x[which.max(dd$y)]
        }
    x.half = x[which (x<=x.mean)]
    null = c (x.half, 2*x.mean - x.half)
    x.sd <- (sd (null))
    pval <- pnorm(x, x.mean, x.sd, lower.tail = FALSE)
    return (pval)
}
