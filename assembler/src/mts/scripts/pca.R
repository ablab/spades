load_clusters <- function(canopy_in, canopy_out, int_contigs) {
  data <- read.table(canopy_in)
  names(data) <- c('contig', sapply(seq(1, dim(data)[2]-1, 1),
                             function(x) {paste('mlt', x, sep='')}))
  binned <- read.table(canopy_out)
  names(binned) <- c('clust', 'contig')
  interesting <- read.table(int_contigs)
  names(interesting) <- c('contig')
  contigs <- merge(x=binned, y=interesting, by='contig')
  droplevels(merge(x=data, y=contigs, by='contig'))
}

do_prc <- function(clusters) {
  prcomp(~ ., data = clusters[, grep('mlt', colnames(clusters))])
}

print_clusters <- function(pr, clust) {
  #depng(filename='tmp.png')
  lev <- levels(factor(clust))
  cols <- 1:length(clust)
  plot(pr$x, col = as.numeric(clust))
  a <- split(as.data.frame(pr$x), clust)
  for (l in lev) {
    x <- a[[l]]
    text(median(x$PC1), median(x$PC2), l)
  }
  legend(0,0, lev, col=cols, pch=1)
  #dev.off()
}

local_data <- function() {
  clusters <- load_clusters("/Volumes/Chihua-Sid/mts/out/sample9.in",
                            "/Volumes/Chihua-Sid/mts/out/sample9.out",
                            "/Volumes/Chihua-Sid/mts/out/70p_3.log")
  
  prc_data <- do_prc(clusters)
  print_clusters(prc_data, clusters$clust)
  prc_data
}
