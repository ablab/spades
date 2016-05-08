combine_pieces <- function(table) {
  res <- table
  res$contig <- sub("_.*", "", res$contig)
  unique(res)
}

load_clusters <- function(canopy_in, canopy_out, int_contigs) {
  data <- read.table(canopy_in)
  names(data) <- c('contig', sapply(seq(1, dim(data)[2]-1, 1),
                             function(x) {paste('mlt', x, sep='')}))
  data <- combine_pieces(data)
  binned <- read.table(canopy_out)
  names(binned) <- c('clust', 'contig')
  binned <- combine_pieces(binned)
  interesting <- read.table(int_contigs)
  names(interesting) <- c('contig', 'alignment')
  contigs <- merge(x=binned, y=interesting, by='contig')
  droplevels(merge(x=data, y=contigs, by='contig'))
}

do_prc <- function(clusters) {
  prcomp(~ ., data = clusters[, grep('mlt', colnames(clusters))])
}

print_clusters <- function(pr, clust, image) {
  if (!missing(image))
    png(filename=image)
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

#For debugging
local_data <- function() {
  clusters <- load_clusters("/Volumes/Chihua-Sid/mts/out/sample9.in",
                            "/Volumes/Chihua-Sid/mts/out/sample9.out",
                            "/Volumes/Chihua-Sid/mts/out/70p_3.log")

  prc_data <- do_prc(clusters)
  print_clusters(prc_data, clusters$clust)
  prc_data
}

args <- commandArgs(trailingOnly = TRUE)
in_fn <- args[1]
out_fn <- args[2]
cont_fn <- args[3]
image_out <- args[4]

clusters <- load_clusters(in_fn, out_fn, cont_fn)
prc_data <- do_prc(clusters)
print_clusters(prc_data, clusters$clust, image_out)
