[![TeamCity Simple Build Status](http://chihua.cab.spbu.ru:3000/app/rest/builds/buildType:(id:SPAdesBasicTests_RunAll)/statusIcon)](http://chihua.cab.spbu.ru:3000/buildConfiguration/SPAdesBasicTests_RunAll?mode=builds)
[![License](https://img.shields.io/badge/licence-GPLv2-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ablab/spades)](https://github.com/ablab/spades/releases/)
[![GitHub Downloads](https://img.shields.io/github/downloads/ablab/spades/total.svg?style=social&logo=github&label=Download)](https://github.com/ablab/spades/releases)
[![BioConda Downloads](https://anaconda.org/bioconda/spades/badges/downloads.svg)](https://anaconda.org/bioconda/spades)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md)


<font size=20>__SPAdes 3.15.5 Manual__</font>

1. [About SPAdes](#sec1) </br>
    1.1. [Supported data types](#sec1.1)</br>
    1.2. [SPAdes pipeline](#sec1.2)</br>
    1.3. [SPAdes performance](#sec1.3)</br>
2. [Installation](#sec2)</br>
    2.1. [Downloading SPAdes Linux binaries](#sec2.1)</br>
    2.2. [Downloading SPAdes binaries for Mac](#sec2.2)</br>
    2.3. [Downloading and compiling SPAdes source code](#sec2.3)</br>
    2.4. [Verifying your installation](#sec2.4)</br>
3. [Running SPAdes](#sec3)</br>
    3.1. [SPAdes input](#sec3.1)</br>
    3.2. [SPAdes command line options](#sec3.2)</br>
    3.3. [Assembling IonTorrent reads](#sec3.3)</br>
    3.4. [Assembling long Illumina paired reads (2x150 and 2x250)](#sec3.4)</br>
    3.5. [HMM-guided mode](#hmm)</br>
    3.6. [SPAdes output](#spadesoutput)</br>
    3.7. [plasmidSPAdes output](#plasmidout)</br>
    3.8. [metaplasmidSPAdes and metaviralSPAdes output](#metapv)</br>
    3.9. [biosyntheticSPAdes output](#bgc)</br>
    3.10. [Assembly evaluation](#eval)</br>
4. [Stand-alone binaries released within SPAdes package](#sec4)</br>
    4.1. [k-mer counting](#sec4.1)</br>
    4.2. [k-mer coverage read filter](#sec4.2)</br>
    4.3. [k-mer cardinality estimating](#sec4.3)</br>
    4.4. [Graph construction](#sec4.4)</br>
    4.5. [Long read to graph alignment](#sec4.5)</br>
        4.5.1. [hybridSPAdes aligner](#sec4.5.1)</br>
        4.5.2. [SPAligner](#sec4.5.2)</br>
5. [Citation](#sec5)</br>
6. [Feedback and bug reports](#sec6)</br>

<a name="sec1"></a>
# About SPAdes

SPAdes - St. Petersburg genome assembler - is an assembly toolkit containing various assembly pipelines. This manual will help you to install and run SPAdes. SPAdes version 3.15.5 was released under GPLv2 on July 14th, 2022 and can be downloaded from <http://cab.spbu.ru/software/spades/>. 

The latest SPAdes paper describing various pipelines in a protocol format is available [here](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102).

<a name="sec1.1"></a>
## Supported data types

The current version of SPAdes works with Illumina or IonTorrent reads and is capable of providing hybrid assemblies using PacBio, Oxford Nanopore and Sanger reads. You can also provide additional contigs that will be used as long reads.

Version 3.15.5 of SPAdes supports paired-end reads, mate-pairs and unpaired reads. SPAdes can take as input several paired-end and mate-pair libraries simultaneously. Note, that SPAdes was initially designed for small genomes. It was tested on bacterial (both single-cell MDA and standard isolates), fungal and other small genomes. SPAdes is not intended for larger genomes (e.g. mammalian size genomes). For such purposes you can use it at your own risk.

If you have high-coverage data for bacterial/viral isolate or multi-cell organism, we highly recommend to use [`--isolate`](#isolate) option.

SPAdes 3.15.5 includes the following additional pipelines:

-   metaSPAdes - a pipeline for metagenomic data sets (see [metaSPAdes options](#meta)).
-   plasmidSPAdes - a pipeline for extracting and assembling plasmids from WGS data sets (see [plasmid options](#plasmid)).
-   metaplasmidSPAdes - a pipeline for extracting and assembling plasmids from *metagenomic* data sets (see [plasmid options](#plasmid)).
-   rnaSPAdes - a *de novo* transcriptome assembler from RNA-Seq data (see [rnaSPAdes manual](assembler/rnaspades_manual.html)).
-   biosyntheticSPAdes - a module for biosynthetic gene cluster assembly with paired-end reads (see [biosynthicSPAdes options](#biosynthetic)).
-   rnaviralSPAdes - a *de novo* assembler tailored for RNA viral datasets (transcriptome, metatranscriptome and metavirome). 
-   coronaSPAdes is a special mode of rnaviralSPAdes specifically aimed for SARS-CoV-2 *de novo* assembly.
-   truSPAdes - (DEPRECATED) a module for TruSeq barcode assembly (see [truSPAdes manual](assembler/truspades_manual.html)).

In addition, we provide several stand-alone binaries with relatively simple command-line interface: [k-mer counting](#sec4.1) (`spades-kmercounter`), [assembly graph construction](#sec4.2) (`spades-gbuilder`) and [long read to graph aligner](#sec4.3) (`spades-gmapper`). To learn options of these tools you can either run them without any parameters or read [this section](#sec4).

[]()

<a name="sec1.2"></a>
## SPAdes pipeline

SPAdes comes in several separate modules:

-   [BayesHammer](http://bioinf.spbau.ru/en/spades/bayeshammer) - read error correction tool for Illumina reads, which works well on both single-cell and standard data sets.
-   IonHammer - read error correction tool for IonTorrent data, which also works on both types of data.
-   SPAdes - iterative short-read genome assembly module; values of K are selected automatically based on the read length and data set type.
-   MismatchCorrector - a tool which improves mismatch and short indel rates in resulting contigs and scaffolds; this module uses the [BWA](http://bio-bwa.sourceforge.net) tool \[[Li H. and Durbin R., 2009](http://www.ncbi.nlm.nih.gov/pubmed/19451168)\]; MismatchCorrector is turned off by default, but we recommend to turn it on (see [SPAdes options section](#correctoropt)).

We recommend to run SPAdes with BayesHammer/IonHammer to obtain high-quality assemblies. However, if you use your own read correction tool, it is possible to turn error correction module off. It is also possible to use only the read error correction stage, if you wish to use another assembler. See the [SPAdes options section](#pipelineopt). []()

<a name="sec1.3"></a>
## SPAdes performance

In this section we give approximate data about SPAdes performance on two data sets:

-   [Standard isolate *E. coli*](https://www.ncbi.nlm.nih.gov/sra/?term=ERR008613); 6.2Gb, 28M reads, 2x100bp, insert size ~ 215bp
-   [MDA single-cell *E. coli*](http://cab.spbu.ru/files/spades_test_datasets/ecoli_sc/); 6.3 Gb, 29M reads, 2x100bp, insert size ~ 270bp (originally downloaded from [here](http://bix.ucsd.edu/projects/singlecell/nbt_data.html))

We ran SPAdes with default parameters using 16 threads on a server with Intel Xeon 2.27GHz processors. BayesHammer runs in approximately half an hour and takes up to 8Gb of RAM to perform read error correction on each data set. Assembly takes about 10 minutes for the *E. coli* isolate data set and 20 minutes for the *E. coli* single-cell data set. Both data sets require about 8Gb of RAM (see notes below). MismatchCorrector runs for about 15 minutes on both data sets, and requires less than 2Gb of RAM. All modules also require additional disk space for storing results (corrected reads, contigs, etc) and temporary files. See the table below for more precise values.

<table border="1" cellpadding="4" cellspacing="0">
<tr>
<td align="right"> Data set </td>
<td colspan="3" align="center"> <i>E. coli</i> isolate </td> 
<td colspan="3" align="center"> <i>E. coli</i> single-cell </td>
</tr>

<tr>
<td> Stage </td>
<td align="center" width="110">  Time  </td> 
<td align="center" width="110">  Peak RAM <br> usage (Gb)  </td>
<td align="center" width="110">  Additional <br> disk space (Gb)  </td>
<td align="center" width="110">  Time </td> 
<td align="center" width="110">  Peak RAM <br> usage (Gb)  </td>
<td align="center" width="110">  Additional <br> disk space (Gb) </td>
</tr>

<tr>
<td> BayesHammer </td>
<td align="center"> 24m </td>
<td align="center"> 7.8 </td>
<td align="center"> 8.5 </td>
<td align="center"> 25m </td>
<td align="center"> 7.7 </td>
<td align="center"> 8.6 </td>
</tr>

<tr>
<td> SPAdes </td>
<td align="center"> 8m </td>
<td align="center"> 8.4 </td>
<td align="center"> 1.4 </td>
<td align="center"> 10m </td>
<td align="center"> 8.3 </td>
<td align="center"> 2.1 </td>
</tr>

<tr>
<td> MismatchCorrector </td>
<td align="center"> 10m </td>
<td align="center"> 1.7 </td>
<td align="center"> 21.4 </td>
<td align="center"> 12m </td>
<td align="center"> 1.8 </td>
<td align="center"> 22.4 </td>
</tr>

<tr>
<td> Whole pipeline </td>
<td align="center"> 42m </td>
<td align="center"> 8.4 </td>
<td align="center"> 23.9 </td>
<td align="center"> 47m </td>
<td align="center"> 8.3 </td>
<td align="center"> 25.1 </td>
</tr>
</table>

Notes:

-   Running SPAdes without preliminary read error correction (e.g. without BayesHammer or IonHammer) will likely require more time and memory.
-   Each module removes its temporary files as soon as it finishes.
-   SPAdes uses 512 Mb per thread for buffers, which results in higher memory consumption. If you set memory limit manually, SPAdes will use smaller buffers and thus less RAM.
-   Performance statistics is given for SPAdes version 3.14.1.

<a name="sec2"></a>
# Installation


SPAdes requires a 64-bit Linux system or Mac OS and Python (supported versions are Python 2.7, and Python3: 3.2 and higher) to be pre-installed on it. To obtain SPAdes you can either download binaries or download source code and compile it yourself. []()

In case of successful installation the following files will be placed in the `bin` directory:

-   `spades.py` (main executable script)
-   `metaspades.py` (main executable script for [metaSPAdes](#meta))
-   `plasmidspades.py` (main executable script for [plasmidSPAdes](#plasmid))
-   `metaplasmidspades.py` (main executable script for [metaplasmidSPAdes](#metaextrachromosomal))
-   `metaviralspades.py` (main executable script for [metaviralSPAdes](#metaextrachromosomal))
-   `rnaspades.py` (main executable script for [rnaSPAdes](rnaspades_manual.html))
-   `truspades.py` (main executable script for [truSPAdes](truspades_manual.html), DEPRECATED)
-   `rnaviralspades.py` (main executable script for rnaviralSPAdes)
-   `coronaspades.py` (wrapper script for coronaSPAdes mode)
-   `spades-core`  (assembly module)
-   `spades-gbuilder`  (standalone graph builder application)
-   `spades-gmapper`  (standalone long read to graph aligner)
-   `spades-kmercount`  (standalone k-mer counting application)
-   `spades-hammer`  (read error correcting module for Illumina reads)
-   `spades-ionhammer`  (read error correcting module for IonTorrent reads)
-   `spades-bwa`  ([BWA](http://bio-bwa.sourceforge.net) alignment module which is required for mismatch correction)
-   `spades-corrector-core`  (mismatch correction module)
-   `spades-truseq-scfcorrection`  (executable used in truSPAdes pipeline)

<a name="sec2.1"></a>
## Downloading SPAdes Linux binaries

To download [SPAdes Linux binaries](http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz) and extract them, go to the directory in which you wish SPAdes to be installed and run:

``` bash

    wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Linux.tar.gz
    tar -xzf SPAdes-3.15.5-Linux.tar.gz
    cd SPAdes-3.15.5-Linux/bin/
```

In this case you do not need to run any installation scripts - SPAdes is ready to use. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

Note, that pre-build binaries do not work on new Linux kernels.

<a name="sec2.2"></a>
## Downloading SPAdes binaries for Mac

To obtain [SPAdes binaries for Mac](http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Darwin.tar.gz), go to the directory in which you wish SPAdes to be installed and run:

``` bash

    curl http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5-Darwin.tar.gz -o SPAdes-3.15.5-Darwin.tar.gz
    tar -zxf SPAdes-3.15.5-Darwin.tar.gz
    cd SPAdes-3.15.5-Darwin/bin/
```

Just as in Linux, SPAdes is ready to use and no further installation steps are required. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

<a name="sec2.3"></a>
## Downloading and compiling SPAdes source code

If you wish to compile SPAdes by yourself you will need the following libraries to be pre-installed:

-   g++ (version 5.3.1 or higher)
-   cmake (version 3.5 or higher)
-   zlib
-   libbz2

If you meet these requirements, you can download the [SPAdes source code](http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5.tar.gz):

``` bash

    wget http://cab.spbu.ru/files/release3.15.5/SPAdes-3.15.5.tar.gz
    tar -xzf SPAdes-3.15.5.tar.gz
    cd SPAdes-3.15.5
```

and build it with the following script:

``` bash

    ./spades_compile.sh
```

SPAdes will be built in the directory `./bin`. If you wish to install SPAdes into another directory, you can specify full path of destination folder by running the following command in `bash` or `sh`:

``` bash

    PREFIX=<destination_dir> ./spades_compile.sh
```

for example:

``` bash

    PREFIX=/usr/local ./spades_compile.sh
```

which will install SPAdes into `/usr/local/bin`.

After installation you will get the same files (listed above) in `./bin` directory (or `<destination_dir>/bin` if you specified PREFIX). We also suggest adding `bin` directory to the `PATH` variable. []()

<a name="sec2.4"></a>
## Verifying your installation

For testing purposes, SPAdes comes with a toy data set (reads that align to first 1000 bp of *E. coli*). To try SPAdes on this data set, run:

``` bash

    <spades installation dir>/bin/spades.py --test
```

If you added `bin` folder from SPAdes installation directory to the `PATH` variable, you can run:

``` bash

    spades.py --test
```

For the simplicity we further assume that `bin` folder from SPAdes installation directory is added to the `PATH` variable.

If the installation is successful, you will find the following information at the end of the log:

``` plain

===== Assembling finished. Used k-mer sizes: 21, 33, 55

 * Corrected reads are in spades_test/corrected/
 * Assembled contigs are in spades_test/contigs.fasta
 * Assembled scaffolds are in spades_test/scaffolds.fasta
 * Assembly graph is in spades_test/assembly_graph.fastg
 * Assembly graph in GFA format is in spades_test/assembly_graph_with_scaffolds.gfa
 * Paths in the assembly graph corresponding to the contigs are in spades_test/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in spades_test/scaffolds.paths

======= SPAdes pipeline finished.

========= TEST PASSED CORRECTLY.

SPAdes log can be found here: spades_test/spades.log

Thank you for using SPAdes!
```

<a name="sec3"></a>
# Running SPAdes

<a name="sec3.1"></a>
<a name="input"></a>
## SPAdes input

SPAdes takes as input paired-end reads, mate-pairs and single (unpaired) reads in FASTA and FASTQ. For IonTorrent data SPAdes also supports unpaired reads in unmapped BAM format (like the one produced by Torrent Server). However, in order to run read error correction, reads should be in FASTQ or BAM format. Sanger, Oxford Nanopore and PacBio CLR reads can be provided in both formats since SPAdes does not run error correction for these types of data.

To run SPAdes 3.15.5 you need at least one library of the following types:

-   Illumina paired-end/high-quality mate-pairs/unpaired reads
-   IonTorrent paired-end/high-quality mate-pairs/unpaired reads
-   PacBio CCS reads

Illumina and IonTorrent libraries should not be assembled together. All other types of input data are compatible. SPAdes should not be used if only PacBio CLR, Oxford Nanopore, Sanger reads or additional contigs are available.

SPAdes supports mate-pair only assembly. However, we recommend to use only high-quality mate-pair libraries in this case (e.g. that do not have a paired-end part). We tested mate-pair only pipeline using Illumina Nextera mate-pairs. See more [here](#hqmp).

Notes:

-   It is strongly suggested to provide multiple paired-end and mate-pair libraries according to their insert size (from smallest to longest).
-   It is not recommended to run SPAdes on PacBio reads with low coverage (less than 5).
-   We suggest not to run SPAdes on PacBio reads for large genomes.
-   SPAdes accepts gzip-compressed files.

<a name="input:pairedend"></a>
### Read-pair libraries

By using command line interface, you can specify up to nine different paired-end libraries, up to nine mate-pair libraries and also up to nine high-quality mate-pair ones. If you wish to use more, you can use [YAML data set file](#inputdata:yaml). We further refer to paired-end and mate-pair libraries simply as to read-pair libraries.

By default, SPAdes assumes that paired-end and high-quality mate-pair reads have forward-reverse (fr) orientation and usual mate-pairs have reverse-forward (rf) orientation. However, different orientations can be set for any library by using SPAdes options.

To distinguish reads in pairs we refer to them as left and right reads. For forward-reverse orientation, the forward reads correspond to the left reads and the reverse reads, to the right. Similarly, in reverse-forward orientation left and right reads correspond to reverse and forward reads, respectively, etc.

Each read-pair library can be stored in several files or several pairs of files. Paired reads can be organized in two different ways:

-   In file pairs. In this case left and right reads are placed in different files and go in the same order in respective files.
-   In interleaved files. In this case, the reads are interlaced, so that each right read goes after the corresponding paired left read.

For example, Illumina produces paired-end reads in two files: `R1.fastq` and `R2.fastq`. If you choose to store reads in file pairs make sure that for every read from `R1.fastq` the corresponding paired read from `R2.fastq` is placed in the respective paired file on the same line number. If you choose to use interleaved files, every read from `R1.fastq` should be followed by the corresponding paired read from `R2.fastq`.

If adapter and/or quality trimming software has been used prior to assembly, files with the orphan reads can be provided as "single read files" for the corresponding read-pair library.

<a name="input:merged"></a>
If you have merged some of the reads from your paired-end (not mate-pair or high-quality mate-pair) library (using tools s.a. [BBMerge](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or [STORM](https://bitbucket.org/yaoornl/align_test/overview)), you should provide the file with resulting reads as a "merged read file" for the corresponding library.
Note that non-empty files with the remaining unmerged left/right reads (separate or interlaced) **must** be provided for the same library (for SPAdes to correctly detect the original read length).

In an unlikely case some of the reads from your mate-pair (or high-quality mate-pair) library are "merged", you should provide the resulting reads as a SEPARATE single-read library.

<a name="input:single"></a>
### Unpaired (single-read) libraries

By using command line interface, you can specify up to nine different single-read libraries. To input more libraries, you can use [YAML data set file](#inputdata:yaml).

Single librairies are assumed to have high quality and a reasonable coverage. For example, you can provide PacBio CCS reads as a single-read library.

Note, that you should not specify PacBio CLR, Sanger reads or additional contigs as single-read libraries, each of them has a separate [option](#inputdata). []()

<a name="input:longreads"></a>
### PacBio and Oxford Nanopore reads

SPAdes can take as an input an unlimited number of PacBio and Oxford Nanopore libraries.

PacBio CLR and Oxford Nanopore reads are used for hybrid assemblies (e.g. with Illumina or IonTorrent). There is no need to pre-correct this kind of data. SPAdes will use PacBio CLR and Oxford Nanopore reads for gap closure and repeat resolution.

For PacBio you just need to have filtered subreads in FASTQ/FASTA format. Provide these filtered subreads using `--pacbio` option. Oxford Nanopore reads are provided with `--nanopore` option.

PacBio CCS/Reads of Insert reads or pre-corrected (using third-party software) PacBio CLR / Oxford Nanopore reads can be simply provided as single reads to SPAdes.

<a name="input:contigs"></a>
### Additional contigs

In case you have contigs of the same genome generated by other assembler(s) and you wish to merge them into SPAdes assembly, you can specify additional contigs using `--trusted-contigs` or `--untrusted-contigs`. First option is used when high quality contigs are available. These contigs will be used for graph construction, gap closure and repeat resolution. Second option is used for less reliable contigs that may have more errors or contigs of unknown quality. These contigs will be used only for gap closure and repeat resolution. The number of additional contigs is unlimited.

Note, that SPAdes does not perform assembly using genomes of closely-related species. Only contigs of the same genome should be specified.

[]()
<a name="sec3.2"></a>
## SPAdes command line options

To run SPAdes from the command line, type

``` bash

    spades.py [options] -o <output_dir>
```

Note that we assume that `bin` forder from SPAdes installation directory is added to the `PATH` variable (provide full path to SPAdes executable otherwise: `<spades installation dir>/bin/spades.py`). []()

<a name="basicopt"></a>
### Basic options

`-o <output_dir> `
    Specify the output directory. Required option.

[]()

<a name="isolate"></a>
`--isolate ` 
    This flag is highly recommended for high-coverage isolate and multi-cell Illumina data; improves the assembly quality and running time. 
    We also recommend to trim your reads prior to the assembly. More details can be found [here](http://cab.spbu.ru/benchmarking-tools-for-de-novo-microbial-assembly/).
    This option is not compatible with `--only-error-correction` or `--careful` options. 

<a name="sc"></a>
`--sc `
    This flag is required for MDA (single-cell) data.

[]()

<a name="meta"></a>
`--meta `   (same as `metaspades.py`)
    This flag is recommended when assembling metagenomic data sets (runs metaSPAdes, see [paper](https://genome.cshlp.org/content/27/5/824.short) for more details). Currently metaSPAdes supports only a **_single_** short-read library which has to be **_paired-end_** (we hope to remove this restriction soon). In addition, you can provide long reads (e.g. using `--pacbio` or `--nanopore` options), but hybrid assembly for metagenomes remains an experimental pipeline and optimal performance is not guaranteed. It does not support [careful mode](#correctoropt) (mismatch correction is not available). In addition, you cannot specify coverage cutoff for metaSPAdes. Note that metaSPAdes might be very sensitive to presence of the technical sequences remaining in the data (most notably adapter readthroughs), please run quality control and pre-process your data accordingly.

[]()

<a name="plasmid"></a>
`--plasmid `   (same as `plasmidspades.py`)
    This flag is required when assembling only plasmids from WGS data sets (runs plasmidSPAdes, see [paper](https://academic.oup.com/bioinformatics/article/32/22/3380/2525610) for the algorithm details). Note, that plasmidSPAdes is not compatible with [single-cell mode](#sc). Additionally, we do not recommend to run plasmidSPAdes on more than one library. 

For plasmidSPAdes output details see [section 3.6](#plasmidout).

[]()

<a name="metaextrachromosomal"></a> 
`--metaplasmid `   (same as `metaplasmidspades.py` and `--meta` `--plasmid`) and

`--metaviral `   (same as `metaviralspades.py`)
 
These options works specially for extracting extrachromosomal elements from metagenomic assemblies. They run similar pipelines that slightly differ in the simplification step; another difference is that for metaviral mode we output linear putative extrachromosomal contigs and for metaplasmid mode we do not.
See [metaplasmid paper](https://genome.cshlp.org/content/29/6/961.short) and [metaviral paper](https://academic.oup.com/bioinformatics/article-abstract/36/14/4126/5837667) for the algorithms details.

For metaplasmidSPAdes/metaviralSPAdes output details see [section 3.7](#metapv).

[]()

Additionally for plasmidSPAdes, metaplasmidSPAdes and metaviralSPAdes we recommend to additionally verify resulting contigs with [viralVerify tool](https://github.com/ablab/viralVerify).

[]()

<a name="biosynthetic"></a>
`--bio `
    This flag is required when assembling only non-ribosomal and polyketide gene clusters from WGS data sets (runs biosyntheticSPAdes, see [paper](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118?top=1) for the algorithm details). biosyntheticSPAdes is supposed to work on isolate or metagenomic WGS dataset. Note, that biosyntheticSPAdes is not compatible with any other modes. See [section 3.8](#bgc) for biosyntheticSPAdes output details.

[]()

<a name="rna"></a>
`--rna `   (same as `rnaspades.py`)
    This flag should be used when assembling RNA-Seq data sets (runs rnaSPAdes). To learn more, see [rnaSPAdes manual](assembler/rnaspades_manual.html).
    Not compatible with `--only-error-correction` or `--careful` options. 

[]()

<a name="rnaviral"></a>
`--rnaviral`   (same as `rnaviralspades.py`)
    This flag should be used when assembling viral RNA-Seq data sets (runs rnaviralSPAdes).
    Not compatible with `--only-error-correction` or `--careful` options. 

`--iontorrent `
    This flag is required when assembling IonTorrent data. Allows BAM files as input. Carefully read [section 3.3](#sec3.3) before using this option.

`--test`
    Runs SPAdes on the toy data set; see [section 2.4](#sec2.4).

`-h` (or `--help`)
    Prints help.

`-v` (or `--version`)
    Prints SPAdes version.

[]()
<a name="pipelineopt"></a>
### Pipeline options

`--only-error-correction`
    Performs read error correction only.

`--only-assembler`
    Runs assembly module only.

[]()
<a name="correctoropt"></a>
`--careful`
    Tries to reduce the number of mismatches and short indels. Also runs MismatchCorrector - a post processing tool, which uses [BWA](http://bio-bwa.sourceforge.net) tool (comes with SPAdes). This option is recommended only for assembly of small genomes. We strongly recommend not to use it for large and medium-size eukaryotic genomes. Note, that this options is is not supported by metaSPAdes and rnaSPAdes. 

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

`--disable-gzip-output`
    Forces read error correction module not to compress the corrected reads. If this options is not set, corrected reads will be in `*.fastq.gz` format.

[]()

<a name="inputdata"></a>
### Input data

<a name="inputdata:pairedend"></a>
#### Specifying single library (paired-end or single-read)

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

<a name="inputdata:multiple"></a>
#### Specifying multiple libraries 

**_Single-read libraries_**

`--s<#> <file_name> `
    File for single-read library number `<#>` (`<#>` = 1,2,..,9). For example, for the first paired-end library the option is: `--s1 <file_name> `
    Do not use `-s` options for single-read libraries, since it specifies unpaired reads for the first paired-end library.

**_Paired-end libraries_**

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
    Should not be confused with FR and RF strand-specificity for RNA-Seq data (see <a href="assembler/rnaspades_manual.html#sec2.3" target="_blank">rnaSPAdes manual</a>). 

<a name="inputdata:matepairs"></a>
**_Mate-pair libraries_**

`--mp<#>-12 <file_name> `
    File with interlaced reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-1 <file_name> `
    File with left reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-2 <file_name> `
    File with right reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9).

`--mp<#>-<or> `
    Orientation of reads for mate-pair library number `<#>` (`<#>` = 1,2,..,9; `<or>` = "fr","rf","ff").
    The default orientation for mate-pair libraries is reverse-forward (`<-- -->`). For example, to specify forward-forward orientation for the first mate-pair library, you should use the flag: `--mp1-ff `

<a name="hqmp"></a>
**_High-quality mate-pair libraries_** (can be used for mate-pair only assembly)

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

<a name="inputdata:longreads"></a>
**_Specifying data for hybrid assembly_**

`--pacbio <file_name> `
    File with PacBio CLR reads. For PacBio CCS reads use `-s` option. More information on PacBio reads is provided in [section 3.1](#input:longreads).

`--nanopore <file_name> `
    File with Oxford Nanopore reads.

`--sanger <file_name> `
    File with Sanger reads

`--trusted-contigs <file_name> `
    Reliable contigs of the same genome, which are likely to have no misassemblies and small rate of other errors (e.g. mismatches and indels). This option is not intended for contigs of the related species.

`--untrusted-contigs <file_name> `
    Contigs of the same genome, quality of which is average or unknown. Contigs of poor quality can be used but may introduce errors in the assembly. This option is also not intended for contigs of the related species.

**_Other input_**

`--assembly-graph <file_name> `
    File with assembly graph. Could only be used in plasmid, metaplasmid, metaviral and biosynthetic mode. The primary purpose of this option to run these pipelines on already constructed and simplified assembly graph this way skipping a large part of SPAdes pipeline. Original reads the graph was constructed from need to be specified as well. Exact k-mer length (via `-k` option) should be provided. Note that the output would be different as compared to standalone runs of these pipelines as they setup graph simplification options as well.


<a name="inputdata:yaml"></a>
**_Specifying input data with YAML data set file (advanced)_**

An alternative way to specify an input data set for SPAdes is to create a [YAML](http://www.yaml.org/) data set file. By using a YAML file you can provide an unlimited number of paired-end, mate-pair and unpaired libraries. Basically, YAML data set file is a text file, in which input libraries are provided as a comma-separated list in square brackets. Each library is provided in braces as a comma-separated list of attributes. The following attributes are available:

-   orientation ("fr", "rf", "ff")
-   type ("paired-end", "mate-pairs", "hq-mate-pairs", "single", "pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs")
-   interlaced reads (comma-separated list of files with interlaced reads)
-   left reads (comma-separated list of files with left reads)
-   right reads (comma-separated list of files with right reads)
-   single reads (comma-separated list of files with single reads or unpaired reads from paired library)
-   merged reads (comma-separated list of files with [merged reads](#input:merged))

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

[]()

<a name="advancedopt"></a>
### Advanced options

`-t <int>` (or `--threads <int>`)
    Number of threads. The default value is 16.

`-m <int>` (or `--memory <int>`)
    Set memory limit in Gb. SPAdes terminates if it reaches this limit. The default value is 250 Gb. Actual amount of consumed RAM will be below this limit. Make sure this value is correct for the given machine. SPAdes uses the limit value to automatically determine the sizes of various buffers, etc.

`--tmp-dir <dir_name>`
    Set directory for temporary files from read error correction. The default value is `<output_dir>/corrected/tmp`

`-k <int,int,...>`
    Comma-separated list of k-mer sizes to be used (all values must be odd, less than 128 and listed in ascending order). If `--sc` is set the default values are 21,33,55. For multicell data sets K values are automatically selected using maximum read length ([see note for assembling long Illumina paired reads for details](#sec3.4)). To properly select K values for IonTorrent data read [section 3.3](#sec3.3).

`--cov-cutoff <float>`
    Read coverage cutoff value. Must be a positive float value, or "auto", or "off". Default value is "off". When set to "auto" SPAdes automatically computes coverage threshold using conservative strategy. Note, that this option is not supported by metaSPAdes.

`--phred-offset <33 or 64>`
    PHRED quality offset for the input reads, can be either 33 or 64. It will be auto-detected if it is not specified.

`--custom-hmms <file or directory>`
    File or directory with amino acid HMMs for [HMM-guided mode](#hmm).


<a name="examples"></a>
### Examples

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

All options for specifying input data can be mixed if needed, but make sure that files for each library are grouped and files with left and right paired reads are listed in the same order. []()

<a name="sec3.3"></a>
## Assembling IonTorrent reads

Only FASTQ or BAM files are supported as input.

The selection of k-mer length is non-trivial for IonTorrent. If the dataset is more or less conventional (good coverage, not high GC, etc), then use our [recommendation for long reads](#sec3.4) (e.g. assemble using k-mer lengths 21,33,55,77,99,127). However, due to increased error rate some changes of k-mer lengths (e.g. selection of shorter ones) may be required. For example, if you ran SPAdes with k-mer lengths 21,33,55,77 and then decided to assemble the same data set using more iterations and larger values of K, you can run SPAdes once again specifying the same output folder and the following options: `--restart-from k77 -k 21,33,55,77,99,127 --mismatch-correction -o <previous_output_dir>`. Do not forget to copy contigs and scaffolds from the previous run. We are planning to tackle issue of selecting k-mer lengths for IonTorrent reads in next versions.

You may need no error correction for Hi-Q enzyme at all. However, we suggest trying to assemble your data with and without error correction and select the best variant.

For non-trivial datasets (e.g. with high GC, low or uneven coverage) we suggest to enable single-cell mode (setting `--sc` option) and use k-mer lengths of 21,33,55. []()

<a name="sec3.4"></a>
<a name="illumina_long"></a>
## Assembling long Illumina paired reads (2x150 and 2x250)

Recent advances in DNA sequencing technology have led to a rapid increase in read length. Nowadays, it is a common situation to have a data set consisting of 2x150 or 2x250 paired-end reads produced by Illumina MiSeq or HiSeq2500. However, the use of longer reads alone will not automatically improve assembly quality. An assembler that can properly take advantage of them is needed.

SPAdes use of iterative k-mer lengths allows benefiting from the full potential of the long paired-end reads. Currently one has to set the assembler options up manually, but we plan to incorporate automatic calculation of necessary options soon.

Please note that in addition to the read length, the insert length also matters a lot. It is not recommended to sequence a 300bp fragment with a pair of 250bp reads. We suggest using 350-500 bp fragments with 2x150 reads and 550-700 bp fragments with 2x250 reads.

<a name="illumina150"></a>
### Multi-cell data set with read length 2x150 bp

Do not turn off SPAdes error correction (BayesHammer module), which is included in SPAdes default pipeline.

If you have enough coverage (50x+), then you may want to try to set k-mer lengths of 21, 33, 55, 77 (selected by default for reads with length 150bp).

Make sure you run assembler with the `--careful` option to minimize number of mismatches in the final contigs.

We recommend that you check the SPAdes log file at the end of the each iteration to control the average coverage of the contigs.

For reads corrected prior to running the assembler:

``` bash

    spades.py -k 21,33,55,77 --careful --only-assembler <your reads> -o spades_output
```

To correct and assemble the reads:

``` bash

    spades.py -k 21,33,55,77 --careful <your reads> -o spades_output
```

<a name="illumina250"></a>
### Multi-cell data set with read lengths 2x250 bp

Do not turn off SPAdes error correction (BayesHammer module), which is included in SPAdes default pipeline.

By default we suggest to increase k-mer lengths in increments of 22 until the k-mer length reaches 127. The exact length of the k-mer depends on the coverage: k-mer length of 127 corresponds to 50x k-mer coverage and higher. For read length 250bp SPAdes automatically chooses K values equal to 21, 33, 55, 77, 99, 127.

Make sure you run assembler with `--careful` option to minimize number of mismatches in the final contigs.

We recommend you to check the SPAdes log file at the end of the each iteration to control the average coverage of the contigs.

For reads corrected prior to running the assembler:

``` bash

    spades.py -k 21,33,55,77,99,127 --careful --only-assembler <your reads> -o spades_output
```

To correct and assemble the reads:

``` bash

    spades.py -k 21,33,55,77,99,127 --careful <your reads> -o spades_output
```

<a name="illumina_long_sc"></a>
### Single-cell data set with read lengths 2 x 150 or 2 x 250

The default k-mer lengths are recommended. For single-cell data sets SPAdes selects k-mer sizes 21, 33 and 55.

However, it might be tricky to fully utilize the advantages of long reads you have. Consider contacting us for more information and to discuss assembly strategy.
[]()

<a name="hmm"></a>
## HMM-guided mode
The majority of SPAdes assembly modes (normal multicell, single-cell, rnaviral, meta and of course biosynthetic) also supports HMM-guided mode as implemented in biosyntheticSPAdes. The detailed description could be found in [biosyntheticSPAdes paper](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118), but in short: amino acid profile HMMs are aligned to the edges of assembly graph. After this the subgraphs containing the set of matches ("domains") are extracted and all possible paths through the domains that are supported both by paired-end data (via scaffolds) and graph topology are obtained (putative biosynthetic gene clusters).

HMM-guided mode could be enabled via providing a set of HMMs via `--custom-hmms` option. In HMM guided mode the set of contigs and scaffolds (see section [SPAdes output](#spadesoutput) for more information about SPAdes output) is kept intact, however additional [biosyntheticSPAdes output](#bgc) represents the output of HMM-guided assembly.

Note that normal biosyntheticSPAdes mode (via `--bio` option) is a bit different from HMM-guided mode: besides using the special set of profile HMMS representing a family of NRSP/PKS domains also includes a set of assembly graph simplification and processing settings aimed for fuller recovery of biosynthetic gene clusters.

Given an increased interest in coronavirus research we developed a coronavirus assembly mode for SPAdes assembler (also known as coronaSPAdes). It allows to assemble full-length coronaviridae genomes from the transcriptomic and metatranscriptomic data. Algorithmically, coronaSPAdes is an rnaviralSPAdes that uses the set of HMMs from [Pfam SARS-CoV-2 2.0](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam_SARS-CoV-2_2.0/) set as well as additional HMMs as outlined by [(Phan et al, 2019)](https://doi.org/10.1093/ve/vey035). coronaSPAdes could be run via a dedicated `coronaspades.py` script. See [coronaSPAdes preprint](https://www.biorxiv.org/content/10.1101/2020.07.28.224584v1) for more information about rnaviralSPAdes,  coronaSPAdes and HMM-guided mode. Output for any HMM-related mode (--bio, --corona, or --custom-hmms flags) is the same with biosyntheticSPAdes' output.

<a name="spadesoutput"></a>
## SPAdes output

SPAdes stores all output files in `<output_dir> `, which is set by the user.

-   `<output_dir>/corrected/` directory contains reads corrected by BayesHammer in `*.fastq.gz` files; if compression is disabled, reads are stored in uncompressed `*.fastq` files
-   `<output_dir>/scaffolds.fasta` contains resulting scaffolds (recommended for use as resulting sequences)
-   `<output_dir>/contigs.fasta` contains resulting contigs
-   `<output_dir>/assembly_graph_with_scaffolds.gfa` contains SPAdes assembly graph and scaffolds paths in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md)
-   `<output_dir>/assembly_graph.fastg` contains SPAdes assembly graph in [FASTG format](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf)
-   `<output_dir>/contigs.paths` contains paths in the assembly graph corresponding to contigs.fasta (see details below)
-   `<output_dir>/scaffolds.paths` contains paths in the assembly graph corresponding to scaffolds.fasta (see details below)

<a name="spadesoutput:contigs"></a>
### Contigs and scaffolds format

Contigs/scaffolds names in SPAdes output FASTA files have the following format:
`>NODE_3_length_237403_cov_243.207`
Here `3` is the number of the contig/scaffold, `237403` is the sequence length in nucleotides and `243.207` is the k-mer coverage for the last (largest) k value used. Note that the k-mer coverage is always lower than the read (per-base) coverage.

In general, SPAdes uses two techniques for joining contigs into scaffolds. First one relies on read pairs and tries to estimate the size of the gap separating contigs. The second one relies on the assembly graph: e.g. if two contigs are separated by a complex tandem repeat, that cannot be resolved exactly, contigs are joined into scaffold with a fixed gap size of 100 bp. Contigs produced by SPAdes do not contain N symbols.

<a name="spadesoutput:graph"></a>
### Assembly graph formats

To view FASTG and GFA files we recommend to use [Bandage visualization tool](http://rrwick.github.io/Bandage/). Note that sequences stored in `assembly_graph.fastg` correspond to contigs before repeat resolution (edges of the assembly graph). Paths corresponding to contigs after repeat resolution (scaffolding) are stored in `contigs.paths` (`scaffolds.paths`) in the format accepted by Bandage (see [Bandage wiki](https://github.com/rrwick/Bandage/wiki/Graph-paths) for details). The example is given below.

Let the contig with the name `NODE_5_length_100000_cov_215.651` consist of the following edges of the assembly graph:

``` plain
    >EDGE_2_length_33280_cov_199.702
    >EDGE_5_length_84_cov_321.414"
    >EDGE_3_length_111_cov_175.304
    >EDGE_5_length_84_cov_321.414"
    >EDGE_4_length_66661_cov_223.548
```

Then, `contigs.paths` will contain the following record:

``` plain
    NODE_5_length_100000_cov_215.651
    2+,5-,3+,5-,4+
```


Since the current version of Bandage does not accept paths with gaps, paths corresponding contigs/scaffolds jumping over a gap in the assembly graph are split by semicolon at the gap positions. For example, the following record

``` plain
    NODE_3_length_237403_cov_243.207
    21-,17-,15+,17-,16+;
    31+,23-,22+,23-,4-
```

states that `NODE_3_length_237403_cov_243.207` corresponds to the path with 10 edges, but jumps over a gap between edges `EDGE_16_length_21503_cov_482.709` and `EDGE_31_length_140767_cov_220.239`.

<a name="spadesoutput:full"></a>
### Complete list of output files

The full list of `<output_dir>` content is presented below:

- scaffolds.fasta - resulting scaffolds (recommended for use as resulting sequences)
- contigs.fasta - resulting contigs
- assembly_graph.fastg - assembly graph
- contigs.paths - contigs paths in the assembly graph
- scaffolds.paths - scaffolds paths in the assembly graph
- before_rr.fasta - contigs before repeat resolution

- corrected/ - files from read error correction
    - configs/ - configuration files for read error correction
    - corrected.yaml - internal configuration file
    - Output files with corrected reads

- params.txt - information about SPAdes parameters in this run
- spades.log - SPAdes log
- dataset.info - internal configuration file
- input_dataset.yaml - internal YAML data set file
- K<##>/ - directory containing intermediate files from the run with K=<##>. These files should not be used as assembly results; use resulting contigs/scaffolds in files mentioned above.


SPAdes will overwrite these files and directories if they exist in the specified `<output_dir>`. []()

<a name="plasmidout"></a>
## plasmidSPAdes output

plasmidSPAdes and metaplasmidSPAdes output only DNA sequences from putative plasmids. Output file names and formats remain the same as in SPAdes (see [previous](#spadesoutput) section), with the following differences.  

For all plasmidSPAdes' contig names in `contigs.fasta`, `scaffolds.fasta` and `assembly_graph.fastg` we append suffix `_component_X`, where `X` is the id of the putative plasmid, which the contig belongs to. Note that plasmidSPAdes may not be able to separate similar plasmids and thus their contigs may appear with the same id. []()  


<a name="metapv"></a>
## metaplasmidSPAdes/metaviralSPAdes output
The repeat resolution and extrachromosomal element detection in metaplasmidSPAdes/metaviralSPAdes is run independently for different coverage cutoffs values (see [paper](https://genome.cshlp.org/content/29/6/961.short) for details). In order to distinguish contigs with putative plasmids detected at different cutoff levels we extend the contig name in FASTA file with cutoff value used for this particular contig (in format `_cutoff_N`). This is why, in the contrast to regular SPAdes pipeline, there might be a contig with `NODE_1_` prefix for each cutoff with potential plasmids detected. In following example, there were detected two potential viruses using cutoff 0, one virus was detected with cutoff 5 and one with cutoff 10.
Also, we add a suffix that shows the structure of the suspective extrachromosomal element.
For metaplasmid mode we output only circular putative plasmids.
For metaviral mode we also output linear putative viruses and linear viruses with simple repeats ('9'-shaped components in the assembly graph) sequences.

``` plain
>NODE_1_length_40003_cov_13.48_cutoff_0_type_circular
>NODE_2_length_30000_cov_4.20_cutoff_0_type_linear
>NODE_1_length_20000_cov_20.42_cutoff_5_type_circular
>NODE_1_length_10000_cov_198.4_cutoff_10_type_linearrepeat
```


<a name="bgc"></a>
## biosyntheticSPAdes output

biosyntheticSPAdes outputs four files of interest:
- scaffolds.fasta – contains DNA sequences from putative biosynthetic gene clusters (BGC). Since each sample may contain multiple BGCs and biosyntheticSPAdes can output several putative DNA sequences for eash cluster, for each contig name we append suffix <code>_cluster_X_candidate_Y</code>, where X is the id of the BGC and Y is the id of the candidate from the BGC.
- raw_scaffolds.fasta – SPAdes scaffolds generated without domain-graph related algorithms. Very close to regular scaffolds.fasta file.
- hmm_statistics.txt – contains statistics about BGC composition in the sample. First, it outputs number of domain hits in the sample. Then, for each BGC candidate we output domain order with positions on the corresponding DNA sequence from scaffolds.fasta.
- domain_graph.dot – contains domain graph structure, that can be used to assess complexity of the sample and structure of BGCs. For more information about domain graph construction, please refer to the paper.

<a name="eval"></a>
## Assembly evaluation

[QUAST](http://cab.spbu.ru/software/quast/) may be used to generate summary statistics (N50, maximum contig length, GC %, \# genes found in a reference list or with built-in gene finding tools, etc.) for a single assembly. It may also be used to compare statistics for multiple assemblies of the same data set (e.g., SPAdes run with different parameters, or several different assemblers).
[]()


<a name="sec4"></a>
# Stand-alone binaries released within SPAdes package

<a name="sec4.1"></a>
## k-mer counting

To provide input data to SPAdes k-mer counting tool `spades-kmercounter ` you may just specify files in [SPAdes-supported formats](#sec3.1) without any flags (after all options) or provide dataset description file in [YAML format](#inputdata:yaml).

Output: <output_dir>/final_kmers - unordered set of kmers in binary format. Kmers from both forward a
nd reverse-complementary reads are taken into account.

Output format: All kmers are written sequentially without any separators. Each kmer takes the same nu
mber of bits. One kmer of length K takes 2*K bits. Kmers are aligned by 64 bits. For example, one kme
r with length=21 takes 8 bytes, with length=33 takes 16 bytes, and with length=55 takes 16 bytes. Eac
h nucleotide is coded with 2 bits: 00 - A, 01 - C, 10 - G, 11 - T.
                                                   
Example:

        For kmer: AGCTCT
        Memory: 6 bits * 2 = 12, 64 bits (8 bytes)
        Let’s describe bytes:
        data[0] = AGCT -> 11 01 10 00 -> 0xd8                                
        data[1] = CT00 -> 00 00 11 01 -> 0x0d
        data[2] = 0000 -> 00 00 00 00 -> 0x00
        data[3] = 0000 -> 00 00 00 00 -> 0x00
        data[4] = 0000 -> 00 00 00 00 -> 0x00
        data[5] = 0000 -> 00 00 00 00 -> 0x00
        data[6] = 0000 -> 00 00 00 00 -> 0x00
        data[7] = 0000 -> 00 00 00 00 -> 0x00

Synopsis: `spades-kmercount [OPTION...] <input files>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format), input files ignored

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-w, --workdir <dir name> `
    working directory to use (default: current directory)

`-b, --bufsize <int> `
    sorting buffer size in bytes, per thread (default 536870912)

`-h, --help `
    print help message


<a name="sec4.2"></a>
## k-mer coverage read filter

`spades-read-filter` is a tool for filtering reads with median kmer coverage less than threshold.

To provide input data to SPAdes k-mer read filter tool `spades-read-filter ` you should provide dataset description file in [YAML format](#inputdata:yaml).

Synopsis: `spades-read-filter [OPTION...] -d <yaml>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format)

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-o, --outdir <dir> `
    output directory to use (default: current directory)

`-c, --cov <value> `
    median kmer count threshold (read pairs, s.t. kmer count median for BOTH reads LESS OR EQUAL to this value will be ignored)

`-h, --help `
    print help message

<a name="sec4.3"></a>
## k-mer cardinality estimating

`spades-kmer-estimating ` is a tool for estimating approximate number of unique k-mers in the provided reads. Kmers from reverse-complementary reads aren"t taken into account for k-mer cardinality estimating.

To provide input data to SPAdes k-mer cardinality estimating tool `spades-kmer-estimating ` you should provide dataset description file in [YAML format](#inputdata:yaml).

Synopsis: `spades-kmer-estimating [OPTION...] -d <yaml>`

The options are:

`-d, --dataset file <file name> `
    dataset description (in YAML format)

`-k, --kmer <int> `
    k-mer length (default: 21)

`-t, --threads <int> `
    number of threads to use (default: number of CPUs)

`-h, --help `
    print help message

<a name="sec4.4"></a>
## Graph construction
Graph construction tool `spades-gbuilder ` has two mandatory options: dataset description file in [YAML format](#inputdata:yaml) and an output file name.

Synopsis: `spades-gbuilder <dataset description (in YAML)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>] [-b <value>] [-unitigs|-fastg|-gfa|-spades]`

Additional options are:

`-k <int> `
    k-mer length used for construction (must be odd)

`-t <int> `
    number of threads

`-tmp-dir <dir_name>  `
    scratch directory to use

`-b <int> `
    sorting buffer size (per thread, in bytes)

`-unitigs `
    k-mer length used for construction (must be odd)

`-fastg `
    output graph in FASTG format

`-gfa `
    output graph in GFA1 format

`-spades `
    output graph in SPAdes internal format


<a name="sec4.5"></a>
## Long read to graph alignment

<a name="sec4.5.1"></a>
### hybridSPAdes aligner
A tool `spades-gmapper ` gives opportunity to extract long read alignments generated with hybridSPAdes pipeline options. It has three mandatory options: dataset description file in [YAML format](#inputdata:yaml), graph file in GFA format and an output file name.

Synopsis: `spades-gmapper <dataset description (in YAML)> <graph (in GFA)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>]`

Additional options are:

`-k <int> `
    k-mer length that was used for graph construction

`-t <int> `
    number of threads

`-tmpdir <dir_name>  `
    scratch directory to use

While `spades-mapper` is a solution for those who work on hybridSPAdes assembly and want to get exactly its intermediate results, [SPAligner](#sec4.5.2) is an end-product application for sequence-to-graph alignment with tunable parameters and output types.  


<a name="sec4.5.2"></a>
### SPAligner
A tool for fast and accurate alignment of nucleotide sequences to assembly graphs. It takes file with sequences (in fasta/fastq format) and assembly in GFA format and outputs long read to graph alignment in various formats (such as tsv, fasta and [GPA](https://github.com/ocxtal/gpa "GPA-format spec")).

Synopsis: `spaligner assembly/src/projects/spaligner_config.yaml -d <value> -s <value> -g <value> -k <value> [-t <value>] [-o <value>]`

Parameters are:

`-d <type> `
    long reads type: nanopore, pacbio

`-s <filename> `
    file with sequences (in fasta/fastq)

`-g <filename> `
    file with graph (in GFA)

`-k <int> `
    k-mer length that was used for graph construction

`-t <int> `
    number of threads (default: 8)

`-o, --outdir <dir> `
    output directory to use (default: spaligner_result/)

For more information on parameters and options please refer to main SPAligner manual (assembler/src/projects/spaligner/README.md).

Also if you want to align protein sequences please refer to our [pre-release version](https://github.com/ablab/spades/releases/tag/spaligner-paper).


<a name="sec5"></a>
# Citation
If you use SPAdes in your research, please cite [our latest paper](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102).

In case you perform hybrid assembly using  PacBio or Nanopore reads, you may also cite [Antipov et al., 2015](http://bioinformatics.oxfordjournals.org/content/early/2015/11/20/bioinformatics.btv688.short). If you use multiple paired-end and/or mate-pair libraries you may additionally cite papers describing SPAdes repeat resolution algorithms [Prjibelski et al., 2014](http://bioinformatics.oxfordjournals.org/content/30/12/i293.short) and [Vasilinetc et al., 2015](http://bioinformatics.oxfordjournals.org/content/31/20/3262.abstract). 

If you use other pipelines, please cite the following papers:

-   metaSPAdes: [Nurk et al., 2017](https://genome.cshlp.org/content/27/5/824.short).
-   plasmidSPAdes: [Antipov et al., 2016](https://academic.oup.com/bioinformatics/article/32/22/3380/2525610).
-   metaplasmidSPAdes / plasmidVerify: [Antipov et al., 2019](https://genome.cshlp.org/content/29/6/961.short)
-   metaviralSPAdes / viralVerify: [Antipov et al., 2020](https://academic.oup.com/bioinformatics/article-abstract/36/14/4126/5837667)
-   rnaSPAdes: [Bushmanova et al., 2019](https://academic.oup.com/gigascience/article/8/9/giz100/5559527).
-   biosyntheticSPAdes: [Meleshko et al., 2019](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118?top=1).
-   coronaSPAdes paper is currently available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.28.224584v1.abstract).

You may also include older papers [Nurk, Bankevich et al., 2013](http://link.springer.com/chapter/10.1007%2F978-3-642-37195-0_13) or [Bankevich, Nurk et al., 2012](http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021), especially if you assemble single-cell data.

In addition, we would like to list your publications that use our software on our website. Please email the reference, the name of your lab, department and institution to [spades.support@cab.spbu.ru](mailto:spades.support@cab.spbu.ru)

<a name="sec6"></a>
# Feedback and bug reports

Your comments, bug reports, and suggestions are very welcomed. They will help us to further improve SPAdes. If you have any troubles running SPAdes, please send us `params.txt` and `spades.log` from the directory `<output_dir>`.

You can leave your comments and bug reports at [our GitHub repository tracker](https://github.com/ablab/spades/issues) or sent it via e-mail: [spades.support@cab.spbu.ru](mailto:spades.support@cab.spbu.ru)


