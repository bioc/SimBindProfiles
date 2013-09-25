filterRegionsPairwise <- function(X, class.group=1, probes=8, probe.max.spacing=200) {

    c <- X[,3]  
    pick.class = which(c == class.group)
    peaks.raw = list()
    p = c()
    for(i in 1:length(pick.class)){
        inx = pick.class[i]
        if(length(p)==0 || (tail(p,1) == inx-1 && (X[inx,2]-X[inx-1,2]) <= probe.max.spacing)){
        p = c(p, inx)
        } else {
            if(length(p) !=0){
                peaks.raw[[length(peaks.raw)+1]] = p
                p = c(inx)
            }
        }
    }
    if(length(p) !=0) {
        peaks.raw[[length(peaks.raw)+1]] = p
    }
    peaks.filtered = peaks.raw[sapply(peaks.raw, length)>=probes]

    if(length(peaks.filtered) >= 1) {
        peaks.bed = data.frame(chr=rep("",length(peaks.filtered)), start=0, end=0, maxLevel=0, meanLevel=0, nProbes=0, stringsAsFactors=FALSE)
        for(i in 1:length(peaks.filtered)){
            p.start = peaks.filtered[[i]][1]
            p.end = tail(peaks.filtered[[i]],1)
            p.nProbes = length(peaks.filtered[[i]])
            if(class.group == 1) {	
                p.maxLevel = max(X[p.start:p.end, 4] - X[p.start:p.end, 5])
                p.meanLevel = mean(X[p.start:p.end, 4] - X[p.start:p.end, 5])
            }
            if(class.group == 2) {	
                p.maxLevel = max(X[p.start:p.end, 5] - X[p.start:p.end, 4])
                p.meanLevel = mean(X[p.start:p.end, 5] - X[p.start:p.end, 4])
            }
            if(class.group == 3) {	
                p.maxLevel = max((X[p.start:p.end, 4] + X[p.start:p.end, 5])/2)
                p.meanLevel = mean((X[p.start:p.end, 4] + X[p.start:p.end, 5])/2)
            }
            peaks.bed$chr[i] = X[p.start, 1]
            peaks.bed$start[i] = X[p.start, 2]
            peaks.bed$end[i] = X[p.end, 2]
            peaks.bed$maxLevel[i] = p.maxLevel
            peaks.bed$meanLevel[i] = p.meanLevel
            peaks.bed$nProbes[i] = p.nProbes
        }
        return(peaks.bed)
    } else {
        cat("Found no regions in class ", class.group, "\n")
    }
}
