fpos <- scan("false_positive.inf")
#perfect <- read.csv(file="perfect.inf",sep=",",head=FALSE)
#hist(fpos,breaks=100)
filt_fpos = fpos[fpos <= 10000]
#for (i in 1:length(fpos)) {
#	if (fpos[i] > 10000) 
#		fpos[i] <- 10000	
#}
jpeg('fpos_hist.jpg')
hist(filt_fpos, breaks=100)
dev.off()

perf <- scan("perfect.inf")
#for (i in 1:length(perf)) {
#    if (perf[i] > 10000)
#        perf[i] <- 10000
#}
filt_perf = perf[perf <= 10000]
jpeg('perf_hist.jpg')
hist(filt_perf, breaks=100)
dev.off()

 superhist2pdf <- function(x, filename = "multi_histograms.pdf",
 dev = "pdf", title = "Superimposed Histograms", nbreaks ="Sturges") {
 junk = NULL
 grouping = NULL
 for(i in 1:length(x)) {
 junk = c(junk,x[[i]])
 grouping <- c(grouping, rep(i,length(x[[i]]))) }
 grouping <- factor(grouping)
 n.gr <- length(table(grouping))
 xr <- range(junk)
 histL <- tapply(junk, grouping, hist, breaks=nbreaks, plot = FALSE)
 maxC <- max(sapply(lapply(histL, "[[", "counts"), max))
 if(dev == "pdf") { pdf(filename, version = "1.4") } else{}
 if((TC <- transparent.cols <- .Device %in% c("pdf", "png"))) {
 cols <- hcl(h = seq(30, by=360 / n.gr, length = n.gr), l = 80, alpha = 0.5) }
 else {
 h.den <- c(10, 15, 20)
 h.ang <- c(45, 15, -30) }
 if(TC) {
 plot(histL[[1]], xlim = xr, ylim= c(0, maxC), col = cols[1], xlab = "x", main = title) }
 else { plot(histL[[1]], xlim = xr, ylim= c(0, maxC), density = h.den[1], angle = h.ang[1], xlab = "x") }
 if(!transparent.cols) {
 for(j in 2:n.gr) plot(histL[[j]], add = TRUE, density = h.den[j], angle = h.ang[j]) } else {
 for(j in 2:n.gr) plot(histL[[j]], add = TRUE, col = cols[j]) }
 invisible()
 if( dev == "pdf") {
 dev.off() }
 }

l = list(filt_fpos,filt_perf)
superhist2pdf(l, nbreaks=100)
