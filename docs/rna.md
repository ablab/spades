# Assembling RNA-Seq data

## About rnaSPAdes

rnaSPAdes is designed primarily for eukaryotic transcriptome assembly from short reads. It can be also used to assemble prokaryotic transcriptomes and meta-transcriptomes. 
rnaSPAdes take as an input at least one paired-end or single-end library. For hybrid assembly you can use PacBio or Oxford Nanopore reads using the same options as in SPAdes (`--pacbio` and `--nanopore` respectively).

- rnaSPAdes does not support `--careful` and `--cov-cutoff` options.
- rnaSPAdes is not compatible with other pipeline modes such as `--meta`, `--sc`, `--plasmid` etc. If you wish to assemble metatranscriptomic data just run rnaSPAdes as it is.
- By default rnaSPAdes uses 2 k-mer sizes, which are automatically detected using read length (approximately one third and half of the maximal read length). We recommend not to change this parameter because smaller k-mer sizes typically result in multiple chimeric (misassembled) transcripts. In case you have any doubts about your run, do not hesitate to contact us via [GitHub issue tracker](https://github.com/ablab/spades/issues).
- Although rnaSPAdes supports IonTorrent reads, it was not sufficiently tested on such kind of data.

## Assembling multiple RNA-Seq libraries
In case you have sequenced several RNA-Seq libraries using the same protocol from different tissues / conditions, and the goal as to assemble a total transcriptome, we suggest to provide all files as a single library (see [SPAdes input options](running.md#input-data)). Note, that sequencing using the same protocol implies that the resulting reads have the same length, insert size and strand-specificity. Transcript quantification for each sample can be done afterwards by separately mapping reads from each library to the assembled transcripts.

When assembling multiple strand-specific libraries, only the first one will be used to determine strand of each transcript. Thus, we suggest not to mix data with different strand-specificity.


## rnaSPAdes-specific options

### Assembling strand-specific data

rnaSPAdes supports strand-specific RNA-Seq datasets. You can set strand-specific type using the following option:

`--ss <type>`
    Use `<type> = rf` when first read in pair corresponds to reverse gene strand (antisense data, e.g. obtained via dUTP protocol) and `<type> = fr` otherwise (forward). 

Note, that strand-specificity is not related and should not be confused with FR and RF orientation of paired reads. RNA-Seq paired-end reads typically have forward-reverse orientation (--> <--), which is assumed by default and no additional options are needed (see [SPAdes input options](running.md#input-data)).

If the data set is single-end use `--ss rf` option when reads are antisense and `--ss fr` otherwise.

### Hybrid transcriptome assembly

rnaSPAdes now supports conventional --pacbio and --nanopore options (see SPAdes manual). Moreover, in addition to long reads you may also provide a separate file with reads capturing the entire transcript sequences using the following options. Full-length transcripts in such reads can be typically detected using the adapters. Note, that FL reads should be trimmed so that the adapters are excluded.

`--fl-rna <file_name>`
    File with PacBio / Nanopore / contigs that capture full-length transcripts.

## rnaSPAdes output
rnaSPAdes outputs one main FASTA file named transcripts.fasta. The corresponding file with paths in the assembly_graph.fastg is transcripts.paths.

In addition rnaSPAdes outputs transcripts with different level of filtration into <output_dir>/:
- `hard_filtered_transcripts.fasta` - includes only long and reliable transcripts with rather high expression.
- `soft_filtered_transcripts.fasta` - includes short and low-expressed transcipts, likely to contain junk sequences.

We reccomend to use main `transcripts.fasta` file in case you don't have any specific needs for you projects. 

Contigs/scaffolds names in rnaSPAdes output FASTA files have the following format:
`>NODE_97_length_6237_cov_11.9819_g8_i2`

Similarly to SPAdes, 97 is the number of the contig, 6237 is its sequence length in nucleotides and 11.9819 is the k-mer coverage. Note that the k-mer coverage is always lower than the read (per-base) coverage. 
g8_i2 correspond to the gene number 8 and isoform number 2 within this gene. Transcripts with the same gene number are presumably received from same or somewhat similar (e.g. paralogous) genes. Note, that the prediction is based on the presence of shared sequences in the transcripts and is very approximate.

## Assembly evaluation
[rnaQUAST](https://github.com/ablab/rnaquast) may be used for transcriptome assembly quality assessment for model organisms when reference genome and gene database are available. rnaQUAST also includes [BUSCO](https://busco.ezlab.org/) and [GeneMarkS-T](http://topaz.gatech.edu/GeneMark/) tools for _de novo_ evaluation.


