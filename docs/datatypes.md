# Tips on SPAdes parameters

## Assembling Illumina genomic data

### Isolated and multi-cell datasets

When assembling standard eukaryotic and bacterial isolated datasets with decent coverage (50x or higher), we strongly recommend to use `--isolate` option.

SPAdes is capable of detecting optimal k-mer sizes automatically. Thus, if the assembly went smoothly without any errors or warnings, there is nothing to worry about.
For example, for read length 100bp the default k values are 21, 33, 55; for 150bp reads SPAdes uses k-mer sizes 21, 33, 55, 77; and for 250bp reads six iterations are used by default: 21, 33, 55, 77, 99, 127.
We strongly recommend *not* to change `-k` parameter unless you are clearly aware about the effect.

### Single-cell data

The default k-mer lengths are recommended. For single-cell data sets SPAdes always selects k-mer sizes 21, 33 and 55.

Do not hesitate to contact us for more information if you plan to assemble single-cell data with long Illumina reads (250 bp and longer).


## Assembling IonTorrent reads

FASTQ or BAM files are supported as input in IonTorrent mode.

The selection of k-mer sizes might be non-trivial for IonTorrent. If the dataset is more or less conventional (good coverage, moderate or low GC, etc), then you can try using a larger k-mer lengths, e.g. 21, 33, 55, 77, 99, 127.

However, due to increased error rate some changes of k-mer lengths (e.g. selection of shorter ones) may be required. For example, if you ran SPAdes with k-mer lengths 21,33,55,77 and then decided to assemble the same data set using more iterations and larger values of K, you can run SPAdes once again specifying the same output folder and the following options: `--restart-from k77 -k 21,33,55,77,99,127 -o <previous_output_dir>`. Do not forget to copy contigs and scaffolds from the previous run.

You may need no error correction for Hi-Q sequencing kit at all. However, we suggest trying to assemble your data with and without error correction and select the best variant.

For non-trivial datasets (e.g. with high GC, low or uneven coverage) we recommend enabling single-cell mode (setting `--sc` option) and use k-mer lengths of 21, 33, 55.
