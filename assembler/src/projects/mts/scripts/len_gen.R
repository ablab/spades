library(dplyr)
library(tidyr)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)
in_fn <- args[1]
out_fn <- args[2]
d<-read.csv(in_fn, header=FALSE, sep="\t")
colnames(d)<-c("name", "bin")
d <- d %>% mutate(sample=str_extract(d$name, "sample\\d+"))
lens <- as.numeric(str_extract(d$name, "(?<=,)\\d+")) - as.numeric(str_extract(d$name, "\\d+(?=,)"))
na_lens <- is.na(lens)
lens[na_lens] <- str_extract(d$name, "(?<=length_)\\d+")[na_lens]
d <- d %>% mutate(len = as.numeric(lens))
#- as.integer(str_extract(d$name, "\\d+,"))
info <- d %>% group_by(bin, sample) %>% summarize(total_len=sum(len))
dth_largest <- function(x, d) {
  x<-sort(x)
  pos<-max(1, length(x) - d)
  x[pos]
}
bin_info<-info %>% group_by(bin) %>% summarize(third_largest = dth_largest(total_len, 3)) %>% arrange(desc(third_largest))
nrow(bin_info %>% filter(third_largest > 1000000))
write.table(bin_info[bin_info$third_largest > 1000000,], out_fn, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
