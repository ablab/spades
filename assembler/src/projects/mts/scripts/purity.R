args <- commandArgs(trailingOnly = TRUE)
stats_dir <- args[1]
out_fn <- args[2]

gfs <- read.table(file.path(stats_dir, "summary", "TSV", "Genome_fraction_(%).tsv"), sep="\t", header=TRUE, row.names=1)
best_cag <- names(gfs)[apply(gfs, 1, which.max)]
names(best_cag) <- rownames(gfs)

stats_table <- read.table(file.path(stats_dir, "combined_reference", "report.tsv"),
                          sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
length <- as.numeric(stats_table["Total length",])
cag_data <- data.frame(length, row.names=names(gfs))

rownames(cag_data) <- names(gfs)

calc_purity <- function(ref) {
    stats_table <- read.table(file.path(stats_dir, "runs_per_reference", ref, "report.tsv"),
                              sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
    cag <- best_cag[[ref]]
    if(!(cag %in% colnames(stats_table)))
        return(0)
    stats <- suppressWarnings(as.numeric(stats_table[,cag]))
    names(stats) <- rownames(stats_table)
    ref_len <- stats[["Reference length"]]
    gf <- stats[["Genome fraction (%)"]] / 100
    al_len <- ref_len * gf
    tot_len <- cag_data[cag, "length"]
    return(al_len / tot_len)
}

purities <- sapply(names(best_cag), function(ref) {
    return(sprintf("%1.2f%% %s", 100 * calc_purity(ref), best_cag[[ref]]))
})

#purities <- data.frame(names(best_cag), purities)

write.table(purities, out_fn, sep="\t", col.names=FALSE, quote=FALSE)
