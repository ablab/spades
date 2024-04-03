# Tips on SPAdes parameters

## Assembling IonTorrent reads

Only FASTQ or BAM files are supported as input.

The selection of k-mer length is non-trivial for IonTorrent. If the dataset is more or less conventional (good coverage, moderate or low GC, etc), then use our [recommendation for long reads](datatypes.md#assembling-long-illumina-paired-reads) (e.g. assemble using k-mer lengths 21,33,55,77,99,127). However, due to increased error rate some changes of k-mer lengths (e.g. selection of shorter ones) may be required. For example, if you ran SPAdes with k-mer lengths 21,33,55,77 and then decided to assemble the same data set using more iterations and larger values of K, you can run SPAdes once again specifying the same output folder and the following options: `--restart-from k77 -k 21,33,55,77,99,127 --mismatch-correction -o <previous_output_dir>`. Do not forget to copy contigs and scaffolds from the previous run.

You may need no error correction for Hi-Q enzyme at all. However, we suggest trying to assemble your data with and without error correction and select the best variant.

For non-trivial datasets (e.g. with high GC, low or uneven coverage) we suggest enabling single-cell mode (setting `--sc` option) and use k-mer lengths of 21,33,55.

## Assembling long Illumina paired reads

Recent advances in DNA sequencing technology have led to a rapid increase in read length. Nowadays, it is a common situation to have a data set consisting of 2x150 or 2x250 paired-end reads produced by Illumina MiSeq or HiSeq2500. However, the use of longer reads alone will not automatically improve assembly quality. An assembler that can properly take advantage of them is needed.

SPAdes use of iterative k-mer lengths allows benefiting from the full potential of the long paired-end reads. Currently one has to set the assembler options up manually, but we plan to incorporate automatic calculation of necessary options soon.

Please note that in addition to the read length, the insert length also matters a lot. It is not recommended to sequence a 300bp fragment with a pair of 250bp reads. We suggest using 350-500 bp fragments with 2x150 reads and 550-700 bp fragments with 2x250 reads.

### Multi-cell data set with read length 2x150 bp

Do not turn off SPAdes error correction (BayesHammer module), which is included in SPAdes default pipeline.

If you have enough coverage (50x+), then you may want to try to set k-mer lengths of 21, 33, 55, 77 (selected by default for reads with length 150bp).

Make sure you run assembler with the `--careful` option to minimize the number of mismatches in the final contigs.

We recommend that you check the SPAdes log file at the end of the each iteration to control the average coverage of the contigs.

For reads corrected prior to running the assembler:

``` bash

    spades.py -k 21,33,55,77 --careful --only-assembler <your reads> -o spades_output
```

To correct and assemble the reads:

``` bash

    spades.py -k 21,33,55,77 --careful <your reads> -o spades_output
```

### Multi-cell data set with read lengths 2x250 bp

Do not turn off SPAdes error correction (BayesHammer module), which is included in SPAdes default pipeline.

By default we suggest increasing k-mer lengths in increments of 22 until the k-mer length reaches 127. The exact length of the k-mer depends on the coverage: k-mer length of 127 corresponds to 50x k-mer coverage and higher. For read length 250bp SPAdes automatically chooses K values equal to 21, 33, 55, 77, 99, 127.

Make sure you run assembler with `--careful` option to minimize the number of mismatches in the final contigs.

We recommend you to check the SPAdes log file at the end of each iteration to control the average coverage of the contigs.

For reads corrected prior to running the assembler:

``` bash

    spades.py -k 21,33,55,77,99,127 --careful --only-assembler <your reads> -o spades_output
```

To correct and assemble the reads:

``` bash

    spades.py -k 21,33,55,77,99,127 --careful <your reads> -o spades_output
```

### Single-cell data set with read lengths 2 x 150 or 2 x 250

The default k-mer lengths are recommended. For single-cell data sets SPAdes selects k-mer sizes 21, 33 and 55.

However, it might be tricky to fully utilize the advantages of long reads you have. Consider contacting us for more information and to discuss assembly strategy.



