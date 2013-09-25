twoGaussiansNull.BF <- function (x, plot=FALSE, var.equal = FALSE,...)
{
    stopifnot(is.numeric(x)) #, p.adj.method %in% p.adjust.methods)
    thisModel <- ifelse(var.equal, "E", "V")
    xmclust <- Mclust(na.omit(x), G = 2, modelNames = thisModel,
        ...)
    nu <- which.min(xmclust$parameters$mean)
    other = setdiff(1:2, nu)

    x.range = range(na.omit(x))
    x2 = seq(x.range[1], x.range[2], length.out=1000)
	if(plot) {
		plot(x2, dnorm(x2, mean=xmclust$parameters$mean[nu], sd=sqrt(xmclust$parameters$variance$sigmasq[nu])), 
		type="l", xlab="signal", ylab="Density", lwd=2)
		lines(x2, dnorm(x2, mean=xmclust$parameters$mean[other], 
		sd=sqrt(xmclust$parameters$variance$sigmasq[other])), col="blue", lwd=2)
		legend("topright", legend=c("Background", "Signal"), lty=c(1,1), 
		col=c("black", "blue"), lwd=c(2,2))
	}

    xp <- pnorm(x, xmclust$parameters$mean[nu], sqrt(xmclust$parameters$variance$sigmasq[nu]),
        lower.tail = FALSE)
    
    return(xp)
}
