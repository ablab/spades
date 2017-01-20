args<-commandArgs(TRUE)
if (length(args) < 3) {                                                                                                                      
  print("Usage: script.R <profiles.in> <binning.out> <bins.prof>")
  quit(save = "no", status = 239)
}

d<-read.table(args[1])
write.table(data.frame(cag="CAG01", name=d$V1), args[2], sep="\t", row.names = F, col.names = F, quote = F)
writeLines(c("CAG01", apply(d[,-1], 2, median)), args[3], sep = "\t")