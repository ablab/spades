library(stringr)

format_ids <- function(table) {
  table$contig <- paste0(str_extract(table$contig, "\\w+\\d+-"), str_replace(str_extract(table$contig, "NODE_\\d+"), "NODE_", ""))
  unique(table)
}

load_binning <- function(canopy_in, canopy_out) {
  data <- read.table(canopy_in)
  names(data) <- c('contig', sapply(seq(1, dim(data)[2]-1, 1),
                                    function(x) {paste('mlt', x, sep='')}))
  data <- format_ids(data)
  binned <- read.table(canopy_out)
  names(binned) <- c('clust', 'contig')
  binned <- format_ids(binned)
  merge(x=data, y=binned, by='contig')
}

load_clusters <- function(canopy_in, canopy_out, int_contigs) {
  data <- load_binning(canopy_in, canopy_out)
  if (missing(int_contigs)) {
    pieces <- split(data, data$clust)[1:10]
    lims <- lapply(pieces, function(x) head(x, 50))
    do.call(rbind, c(lims, list(make.row.names=FALSE)))
  } else {
    interesting <- read.table(int_contigs)
    names(interesting) <- c('contig', 'length', 'alignment', 'ref')
    droplevels(merge(x=data, y=interesting, by='contig'))
  }
}

do_prc <- function(clusters) {
  prcomp(~ ., data = clusters[, grep('mlt', colnames(clusters))])
}

print_clusters <- function(pr, clust, image) {
  if (!missing(image))
    png(filename=image, width=1024, height=768)
  lev <- levels(factor(clust))
  cols <- 1:length(lev)
  #layout(rbind(1,2), heights=c(7,1))
  plot(pr$x, col = as.numeric(clust), xlim=c(-100, 200), ylim=c(-50,50))
  a <- split(as.data.frame(pr$x), clust)
  for (l in lev) {
    x <- a[[l]]
    text(median(x$PC1), median(x$PC2), l)
  }
  legend("center", "bottom", legend=lev, col=cols, pch=1)
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
if (length(args) < 4) {
  image_out <- args[3]
  clusters <- load_clusters(in_fn, out_fn)
} else {
  cont_fn <- args[3]
  image_out <- args[4]
  clusters <- load_clusters(in_fn, out_fn, cont_fn)
}

print(clusters[1:10,])
prc_data <- do_prc(clusters)
print_clusters(prc_data, clusters$clust, image_out)
