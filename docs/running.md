# SPAdes command line options

To run SPAdes from the command line, type

``` bash

    spades.py [options] -o <output_dir>
```

Note that we assume that the `bin` folder from SPAdes installation directory is added to the `PATH` variable (provide full path to SPAdes executable otherwise: `<spades installation dir>/bin/spades.py`).

## Basic options and modes

`-o <output_dir> `
    Specify the output directory. Required option.

`--isolate `
    This flag is highly recommended for high-coverage isolate and multi-cell Illumina data; improves the assembly quality and running time.
    We also recommend trimming your reads prior to the assembly.
    This option is not compatible with `--only-error-correction` or `--careful` options.

`--sc `
    This flag is required for MDA (single-cell) data.

`--meta `   (same as `metaspades.py`)
    This flag is recommended when assembling metagenomic data sets (runs metaSPAdes, see [paper](https://genome.cshlp.org/content/27/5/824.short) for more details). Currently metaSPAdes supports only a **_single_** short-read library which has to be **_paired-end_** (we hope to remove this restriction soon). In addition, you can provide long reads (e.g. using `--pacbio` or `--nanopore` options), but hybrid assembly for metagenomes remains an experimental pipeline and optimal performance is not guaranteed. It does not support [careful mode](running.md#pipeline-options) (mismatch correction is not available). In addition, you cannot specify coverage cutoff for metaSPAdes. Note that metaSPAdes might be very sensitive to the presence of the technical sequences remaining in the data (most notably adapter readthroughs), please run quality control and pre-process your data accordingly.

`--plasmid `   (same as `plasmidspades.py`)
    This flag is required when assembling only plasmids from WGS data sets (runs plasmidSPAdes, see [paper](https://academic.oup.com/bioinformatics/article/32/22/3380/2525610) for the algorithm details). Note, that plasmidSPAdes is not compatible with single-cell mode (`--sc`). Additionally, we do not recommend to run plasmidSPAdes on more than one library.

See [plasmidSPAdes output section](output.md#plasmidspades-output) for details.

`--metaplasmid `   (same as `metaplasmidspades.py` and `--meta` `--plasmid`) and

`--metaviral `   (same as `metaviralspades.py`)

These options work specially for extracting extrachromosomal elements from metagenomic assemblies. They run similar pipelines that slightly differ in the simplification step; another difference is that for metaviral mode we output linear putative extrachromosomal contigs and for metaplasmid mode we do not.
See [metaplasmid paper](https://genome.cshlp.org/content/29/6/961.short) and [metaviral paper](https://academic.oup.com/bioinformatics/article-abstract/36/14/4126/5837667) for the algorithms details.

See [metaplasmidSPAdes/metaviralSPAdes section](output.md#metaplasmidspades-and-metaviralspades-output) for details see.

Additionally for plasmidSPAdes, metaplasmidSPAdes and metaviralSPAdes we recommend verifying resulting contigs with [viralVerify tool](https://github.com/ablab/viralVerify).

`--bio `
    This flag is required when assembling only non-ribosomal and polyketide gene clusters from WGS data sets (runs biosyntheticSPAdes, see [paper](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118?top=1) for the algorithm details). biosyntheticSPAdes is supposed to work on isolated or metagenomic WGS dataset. Note, that biosyntheticSPAdes is not compatible with any other modes. See [biosyntheticSPAdes output section](output.md#biosyntheticspades-output) for details.

`--rna `   (same as `rnaspades.py`)
    This flag should be used when assembling RNA-Seq data sets (runs rnaSPAdes). To learn more, see [rnaSPAdes manual](rna.md).
    Not compatible with `--only-error-correction` or `--careful` options.

`--rnaviral`   (same as `rnaviralspades.py`)
    This flag should be used when assembling viral RNA-Seq data sets (runs rnaviralSPAdes).
    Not compatible with `--only-error-correction` or `--careful` options.

`--iontorrent `
    This flag is required when assembling IonTorrent data. Allows BAM files as input. Carefully read [IonTorrent section](datatypes.md#assembling-iontorrent-reads) before using this option.

`--test`
    Runs SPAdes on the toy data set; see [installation](installation.md#verifying-your-installation) for details.

`-h` (or `--help`)
    Prints help.

`-v` (or `--version`)
    Prints SPAdes version.


## Pipeline options

`--only-error-correction`
    Performs read error correction only.

`--only-assembler`
    Runs assembly module only.

`--careful`
    Tries to reduce the number of mismatches and short indels. Also runs MismatchCorrector - a post processing tool, which uses [BWA](http://bio-bwa.sourceforge.net) tool (comes with SPAdes). This option is recommended only for assembly of small genomes. We strongly recommend not to use it for large and medium-size eukaryotic genomes. Note that this option is not supported by metaSPAdes and rnaSPAdes.

`--continue`
    Continues SPAdes run from the specified output folder starting from the last available check-point. Check-points are made after:

-   error correction module is finished
-   iteration for each specified K value of assembly module is finished
-   mismatch correction is finished for contigs or scaffolds

For example, if specified K values are 21, 33 and 55 and SPAdes was stopped or crashed during assembly stage with K = 55, you can run SPAdes with the `--continue` option specifying the same output directory. SPAdes will continue the run starting from the assembly stage with K = 55. Error correction module and iterations for K equal to 21 and 33 will not be run again. If `--continue` is set, the only allowed option is `-o <output_dir> `.

`--restart-from <check_point>`
    Restart SPAdes run from the specified output folder starting from the specified check-point. Check-points are:

-   `ec` - start from error correction
-   `as` - restart assembly module from the first iteration
-   `k<int>` - restart from the iteration with specified k values, e.g. k55 (not available in RNA-Seq mode)
-   `mc` - restart mismatch correction
-   `last` - restart from the last available check-point (similar to `--continue`)

In contrast to the `--continue` option, you can change some of the options when using `--restart-from`. You can change any option except: all basic options, all options for specifying input data (including `--dataset`), `--only-error-correction` option and `--only-assembler` option. For example, if you ran assembler with k values 21,33,55 without mismatch correction, you can add one more iteration with k=77 and run mismatch correction step by running SPAdes with following options:
`--restart-from k55 -k 21,33,55,77 --mismatch-correction -o <previous_output_dir>`.
Since all files will be overwritten, do not forget to copy your assembly from the previous run if you need it.

`--checkpoints <mode>`
    Make intermediate checkpoints that will allow restarting SPAdes from an internal stage. Available modes are `none` (default), `all` (makes all possible checkpoints) and `last` (makes a checkpoint only if SPAdes crashes).

Note:
- this option is NOT mandatory for using `--restart-from` and `--continue` options, but may speed them up;
- making checkpoints may take more time and a significant amount of disk space.

`--disable-gzip-output`
    Forces read error correction module not to compress the corrected reads. If this options is not set, corrected reads will be in `*.fastq.gz` format.


## Input data

### Specifying single library

`--12 <file_name> `
    File with interlaced forward and reverse paired-end reads.

`-1 <file_name> `
    File with forward reads.

`-2 <file_name> `
    File with reverse reads.

`--merged <file_name> `
    File with merged paired reads.
    If the properties of the library permit, overlapping paired-end reads can be merged using special software.
    Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) **must** be provided for the same library for SPAdes to correctly detect the original read length.

`-s <file_name> `
    File with unpaired reads.

### Specifying multiple libraries

#### Single-read libraries

`--s<#> <file_name> `
    File for single-read library number `<#>` (`<#>` = 1,2,..,9). For example, for the first paired-end library the option is: `--s1 <file_name> `
    Do not use `-s` options for single-read libraries, since it specifies unpaired reads for the first paired-end library.

#### Paired-end libraries

`--pe<#>-12 <file_name> `
    File with interlaced reads for paired-end library number `<#>` (`<#>` = 1,2,..,9). For example, for the first single-read library the option is: `--pe1-12 <file_name> `

`--pe<#>-1 <file_name> `
    File with left reads for paired-end library number `<#>` (`<#>` = 1,2,..,9).

`--pe<#>-2 <file_name> `
    File with right reads for paired-end library number `<#>` (`<#>` = 1,2,..,9).

`--pe<#>-m <file_name> `
    File with merged reads from paired-end library number `<#>` (`<#>` = 1,2,..,9)
    If the properties of the library permit, paired reads can be merged using special software. Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) **must** be provided for the same library for SPAdes to correctly detect the original read length.

`--pe<#>-s <file_name> `
    File with unpaired reads from paired-end library number `<#>` (`<#>` = 1,2,..,9)
    For example, paired reads can become unpaired during the error correction procedure.

`--pe<#>-<or> `
    Orientation of reads for paired-end library number `<#>` (`<#>` = 1,2,..,9; `<or>` = "fr","rf","ff").
    The default orientation for paired-end libraries is forward-reverse (`--> <--`). For example, to specify reverse-forward orientation for the second paired-end library, you should use the flag: `--pe2-rf `
    Should not be confused with FR and RF strand-specificity for RNA-Seq data (see [rnaSPAdes manual](rna.md)).

#### Mate-pair libraries

`--mp<#>-12 <file_name> `
    File with interlaced reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-1 <file_name> `
    File with left reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-2 <file_name> `
    File with right reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-<or> `
    Orientation of reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9; `<or>` = "fr","rf","ff").
    The default orientation for mate-pair libraries is reverse-forward (`<-- -->`). For example, to specify forward-forward orientation for the first mate-pair library, you should use the flag: `--mp1-ff `

#### High-quality mate-pair libraries

High-quality MP data can be used for mate-pair only assembly.

`--hqmp<#>-12 <file_name> `
    File with interlaced reads for high-quality mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--hqmp<#>-1 <file_name> `
    File with left reads for high-quality mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--hqmp<#>-2 <file_name> `
    File with right reads for high-quality mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--hqmp<#>-s <file_name> `
    File with unpaired reads from high-quality mate-pair library number `<#>` (`<#>` = 1,2,..,9)

`--hqmp<#>-<or> `
    Orientation of reads for high-quality mate-pair library number `<#>` (`<#>` = 1,2,..,9; `<or>` = "fr","rf","ff").
    The default orientation for high-quality mate-pair libraries is forward-reverse (`--> <--`). For example, to specify reverse-forward orientation for the first high-quality mate-pair library, you should use the flag: `--hqmp1-rf `

#### Specifying data for hybrid assembly

`--pacbio <file_name> `
    File with PacBio CLR reads. For PacBio CCS reads use `-s` option. More information on PacBio reads is provided [here](hybrid.md#pacbio-and-oxford-nanopore-reads).

`--nanopore <file_name> `
    File with Oxford Nanopore reads.

`--sanger <file_name> `
    File with Sanger reads

`--trusted-contigs <file_name> `
    Reliable contigs of the same genome, which are likely to have no misassemblies and small rate of other errors (e.g. mismatches and indels). This option is not intended for contigs of the related species.

`--untrusted-contigs <file_name> `
    Contigs of the same genome, quality of which is average or unknown. Contigs of poor quality can be used but may introduce errors in the assembly. This option is also not intended for contigs of the related species.

#### Other input

`--assembly-graph <file_name> `
    File with assembly graph. Could only be used in plasmid, metaplasmid, metaviral and biosynthetic mode. The primary purpose of this option is to run these pipelines on already constructed and simplified assembly graphs, thus skipping a large part of SPAdes pipeline. Original reads the graph was constructed from need to be specified as well. Exact k-mer length (via `-k` option) should be provided. Note that the output would be different as compared to standalone runs of these pipelines as they set up graph simplification options as well.


### Specifying multiple libraries with YAML data set file

An alternative way to specify an input data set for SPAdes is to create a [YAML](http://www.yaml.org/) data set file. By using a YAML file you can provide an unlimited number of paired-end, mate-pair and unpaired libraries. Basically, a YAML data set file is a text file, in which input libraries are provided as a comma-separated list in square brackets. Each library is provided in braces as a comma-separated list of attributes. The following attributes are available:

-   orientation ("fr", "rf", "ff")
-   type ("paired-end", "mate-pairs", "hq-mate-pairs", "single", "pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs")
-   interlaced reads (comma-separated list of files with interlaced reads)
-   left reads (comma-separated list of files with left reads)
-   right reads (comma-separated list of files with right reads)
-   single reads (comma-separated list of files with single reads or unpaired reads from paired library)
-   merged reads (comma-separated list of files with [merged reads](input.md#read-pair-libraries))

To properly specify a library you should provide its type and at least one file with reads. Orientation is an optional attribute. Its default value is "fr" (forward-reverse) for paired-end libraries and "rf" (reverse-forward) for mate-pair libraries.

The value for each attribute is given after a colon. Comma-separated lists of files should be given in square brackets. For each file you should provide its full path in double quotes. Make sure that files with right reads are given in the same order as corresponding files with left reads.

For example, if you have one paired-end library split into two pairs of files:

``` bash

    lib_pe1_left_1.fastq
    lib_pe1_right_1.fastq
    lib_pe1_left_2.fastq
    lib_pe1_right_2.fastq
```

one mate-pair library:

``` bash

    lib_mp1_left.fastq
    lib_mp1_right.fastq
```

and PacBio CCS and CLR reads:

``` bash

    pacbio_ccs.fastq
    pacbio_clr.fastq
```

YAML file should look like this:

``` bash

    [
      {
        orientation: "fr",
        type: "paired-end",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_pe1_right_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_pe1_right_2.fastq"
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_pe1_left_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_pe1_left_2.fastq"
        ]
      },
      {
        orientation: "rf",
        type: "mate-pairs",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_mp1_right.fastq"
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_mp1_left.fastq"
        ]
      },
      {
        type: "single",
        single reads: [
          "/FULL_PATH_TO_DATASET/pacbio_ccs.fastq"
        ]
      },
      {
        type: "pacbio",
        single reads: [
          "/FULL_PATH_TO_DATASET/pacbio_clr.fastq"
        ]
      }
    ]
```

Once you have created a YAML file save it with `.yaml` extension (e.g. as `my_data_set.yaml`) and run SPAdes using the `--dataset` option:
`--dataset <your YAML file>`
Notes:

-   The `--dataset` option cannot be used with any other options for specifying input data.
-   We recommend to nest all files with long reads of the same data type in a single library block.



## Advanced options

`-t <int>` (or `--threads <int>`)
    Number of threads. The default value is 16.

`-m <int>` (or `--memory <int>`)
    Set memory limit in Gb. SPAdes terminates if it reaches this limit. The default value is 250 Gb. Actual amount of RAM consumed will be below this limit. Make sure this value is correct for the given machine. SPAdes uses the limit value to automatically determine the sizes of various buffers, etc.

`--tmp-dir <dir_name>`
    Set directory for temporary files from read error correction. The default value is `<output_dir>/corrected/tmp`

`-k <int,int,...>`
    Comma-separated list of k-mer sizes to be used (all values must be odd, less than 128 and listed in ascending order). If `--sc` is set the default values are 21,33,55. For multicell data sets K values are automatically selected using maximum read length ([see note for assembling long Illumina paired reads for details](datatypes.md#assembling-long-illumina-paired-reads)). To properly select K values for IonTorrent data read [this section](datatypes.md#assembling-iontorrent-reads).

`--cov-cutoff <float>`
    Read coverage cutoff value. Must be a positive float value, or "auto", or "off". Default value is "off". When set to "auto" SPAdes automatically computes coverage threshold using conservative strategy. Note, that this option is not supported by metaSPAdes.

`--phred-offset <33 or 64>`
    PHRED quality offset for the input reads, can be either 33 or 64. It will be auto-detected if it is not specified.

`--custom-hmms <file or directory>`
    File or directory with amino acid HMMs for [HMM-guided mode](hmm.md#hmm-guided-mode).


## Examples

To test the toy data set, you can also run the following command from the SPAdes `bin` directory:

``` bash

    spades.py --pe1-1 ../share/spades/test_dataset/ecoli_1K_1.fq.gz \
    --pe1-2 ../share/spades/test_dataset/ecoli_1K_2.fq.gz -o spades_test
```

If you have your library separated into several pairs of files, for example:

``` bash

    lib1_forward_1.fastq
    lib1_reverse_1.fastq
    lib1_forward_2.fastq
    lib1_reverse_2.fastq
```

make sure that corresponding files are given in the same order:

``` bash

    spades.py --pe1-1 lib1_forward_1.fastq --pe1-2 lib1_reverse_1.fastq \
    --pe1-1 lib1_forward_2.fastq --pe1-2 lib1_reverse_2.fastq \
    -o spades_output
```

Files with interlacing paired-end reads or files with unpaired reads can be specified in any order with one file per option, for example:

``` bash

    spades.py --pe1-12 lib1_1.fastq --pe1-12 lib1_2.fastq \
    --pe1-s lib1_unpaired_1.fastq --pe1-s lib1_unpaired_2.fastq \
    -o spades_output
```

If you have several paired-end and mate-pair reads, for example:

paired-end library 1

``` bash

    lib_pe1_left.fastq
    lib_pe1_right.fastq
```

mate-pair library 1

``` bash

    lib_mp1_left.fastq
    lib_mp1_right.fastq
```

mate-pair library 2

``` bash

    lib_mp2_left.fastq
    lib_mp2_right.fastq
```

make sure that files corresponding to each library are grouped together:

``` bash

    spades.py --pe1-1 lib_pe1_left.fastq --pe1-2 lib_pe1_right.fastq \
    --mp1-1 lib_mp1_left.fastq --mp1-2 lib_mp1_right.fastq \
    --mp2-1 lib_mp2_left.fastq --mp2-2 lib_mp2_right.fastq \
    -o spades_output
```

If you have IonTorrent unpaired reads, PacBio CLR and additional reliable contigs:

``` bash

    it_reads.fastq
    pacbio_clr.fastq
    contigs.fasta
```

run SPAdes with the following command:

``` bash

    spades.py --iontorrent -s it_reads.fastq \
    --pacbio pacbio_clr.fastq --trusted-contigs contigs.fastq \
    -o spades_output
```

If a single-read library is split into several files:

``` bash

    unpaired1_1.fastq
    unpaired1_2.fastq
    unpaired1_3.fasta
```

specify them as one library:

``` bash

    spades.py --s1 unpaired1_1.fastq \
    --s1 unpaired1_2.fastq --s1 unpaired1_3.fastq \
    -o spades_output
```

All options for specifying input data can be mixed if needed, but make sure that files for each library are grouped and files with left and right paired reads are listed in the same order.

