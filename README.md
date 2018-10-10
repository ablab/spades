<font size=20>__SPAdes 3.13.0 Manual__</font>


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
    3.5. [SPAdes output](#sec3.5)</br>
    3.6. [plasmidSPAdes output](#sec3.6)</br>
    3.7. [Assembly evaluation](#sec3.7)</br>
4. [Stand-alone binaries released within SPAdes package](#sec4)</br>
    4.1. [k-mer counting](#sec4.1)</br>
    4.2. [Graph construction](#sec4.2)</br>
    4.3. [Long read to graph aligner](#sec4.3)</br>
5. [Citation](#sec5)</br>
6. [Feedback and bug reports](#sec6)</br>

<a name="sec1"></a>
# About SPAdes

SPAdes – St. Petersburg genome assembler – is an assembly toolkit containing various assembly pipelines. This manual will help you to install and run SPAdes. SPAdes version 3.13.0 was released under GPLv2 on October 11, 2018 and can be downloaded from <http://cab.spbu.ru/software/spades/>. []()

<a name="sec1.1"></a>
## Supported data types

The current version of SPAdes works with Illumina or IonTorrent reads and is capable of providing hybrid assemblies using PacBio, Oxford Nanopore and Sanger reads. You can also provide additional contigs that will be used as long reads.

Version 3.13.0 of SPAdes supports paired-end reads, mate-pairs and unpaired reads. SPAdes can take as input several paired-end and mate-pair libraries simultaneously. Note, that SPAdes was initially designed for small genomes. It was tested on bacterial (both single-cell MDA and standard isolates), fungal and other small genomes. SPAdes is not intended for larger genomes (e.g. mammalian size genomes). For such purposes you can use it at your own risk.

SPAdes 3.13.0 includes the following additional pipelines:

-   metaSPAdes – a pipeline for metagenomic data sets (see [metaSPAdes options](#meta)).
-   plasmidSPAdes – a pipeline for extracting and assembling plasmids from WGS data sets (see [plasmidSPAdes options](#plasmid)).
-   rnaSPAdes – a *de novo* transcriptome assembler from RNA-Seq data (see [rnaSPAdes manual](rnaspades_manual.html)).
-   truSPAdes – a module for TruSeq barcode assembly (see [truSPAdes manual](truspades_manual.html)).

In addition, we provide several stand-alone binaries with relatively simple command-line interface: [k-mer counting](#sec4.1) (`spades-kmercounter`), [assembly graph construction](#sec4.2) (`spades-gbuilder`) and [long read to graph aligner](#sec4.3) (`spades-gmapper`). To learn options of these tools you can either run them without any parameters or read [this section](#sec4).

[]()

<a name="sec1.2"></a>
## SPAdes pipeline

SPAdes comes in several separate modules:

-   [BayesHammer](http://bioinf.spbau.ru/en/spades/bayeshammer) – read error correction tool for Illumina reads, which works well on both single-cell and standard data sets.
-   IonHammer – read error correction tool for IonTorrent data, which also works on both types of data.
-   SPAdes – iterative short-read genome assembly module; values of K are selected automatically based on the read length and data set type.
-   MismatchCorrector – a tool which improves mismatch and short indel rates in resulting contigs and scaffolds; this module uses the [BWA](http://bio-bwa.sourceforge.net) tool \[[Li H. and Durbin R., 2009](http://www.ncbi.nlm.nih.gov/pubmed/19451168)\]; MismatchCorrector is turned off by default, but we recommend to turn it on (see [SPAdes options section](#correctoropt)).

We recommend to run SPAdes with BayesHammer/IonHammer to obtain high-quality assemblies. However, if you use your own read correction tool, it is possible to turn error correction module off. It is also possible to use only the read error correction stage, if you wish to use another assembler. See the [SPAdes options section](#pipelineopt). []()

<a name="sec1.3"></a>
## SPAdes' performance

In this section we give approximate data about SPAdes' performance on two data sets:

-   [Standard isolate *E. coli*](http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_mc/); 6.2Gb, 28M reads, 2x100bp, insert size ~ 215bp
-   [MDA single-cell *E. coli*](http://spades.bioinf.spbau.ru/spades_test_datasets/ecoli_sc/); 6.3 Gb, 29M reads, 2x100bp, insert size ~ 270bp

We ran SPAdes with default parameters using 16 threads on a server with Intel Xeon 2.27GHz processors. BayesHammer runs in approximately half an hour and takes up to 8Gb of RAM to perform read error correction on each data set. Assembly takes about 10 minutes for the *E. coli* isolate data set and 20 minutes for the *E. coli* single-cell data set. Both data sets require about 8Gb of RAM (see notes below). MismatchCorrector runs for about 15 minutes on both data sets, and requires less than 2Gb of RAM. All modules also require additional disk space for storing results (corrected reads, contigs, etc) and temporary files. See the table below for more precise values.

<table border="1" cellpadding="4" cellspacing="0">
<tr>
<td align="right"> Data set &nbsp; </td>
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
-   Performance statistics is given for SPAdes version 3.13.0.

<a name="sec2"></a>
#Installation


SPAdes requires a 64-bit Linux system or Mac OS and Python (supported versions are Python2: 2.4–2.7, and Python3: 3.2 and higher) to be pre-installed on it. To obtain SPAdes you can either download binaries or download source code and compile it yourself. []()

In case of successful installation the following files will be placed in the `bin` directory:

-   `spades.py` (main executable script)
-   `metaspades.py` (main executable script for [metaSPAdes](#meta))
-   `plasmidspades.py` (main executable script for [plasmidSPAdes](#plasmid))
-   `rnaspades.py` (main executable script for [rnaSPAdes](rnaspades_manual.html))
-   `truspades.py` (main executable script for [truSPAdes](truspades_manual.html))
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

To download [SPAdes Linux binaries](http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz) and extract them, go to the directory in which you wish SPAdes to be installed and run:

``` bash

    wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
    tar -xzf SPAdes-3.13.0-Linux.tar.gz
    cd SPAdes-3.13.0-Linux/bin/
```

In this case you do not need to run any installation scripts – SPAdes is ready to use. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

<a name="sec2.2"></a>
## Downloading SPAdes binaries for Mac

To obtain [SPAdes binaries for Mac](http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Darwin.tar.gz), go to the directory in which you wish SPAdes to be installed and run:

``` bash

    curl http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Darwin.tar.gz -o SPAdes-3.13.0-Darwin.tar.gz
    tar -zxf SPAdes-3.13.0-Darwin.tar.gz
    cd SPAdes-3.13.0-Darwin/bin/
```

Just as in Linux, SPAdes is ready to use and no further installation steps are required. We also suggest adding SPAdes installation directory to the `PATH` variable. []()

<a name="sec2.3"></a>
## Downloading and compiling SPAdes source code

If you wish to compile SPAdes by yourself you will need the following libraries to be pre-installed:

-   g++ (version 5.3.1 or higher)
-   cmake (version 2.8.12 or higher)
-   zlib
-   libbz2

If you meet these requirements, you can download the [SPAdes source code](http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz):

``` bash

    wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0.tar.gz
    tar -xzf SPAdes-3.13.0.tar.gz
    cd SPAdes-3.13.0
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

After installation you will get the same files (listed above) in `./bin` directory (or `<destination_dir>/bin` if you specified PREFIX). We also suggest adding SPAdes installation directory to the `PATH` variable. []()

<a name="sec2.4"></a>
## Verifying your installation

For testing purposes, SPAdes comes with a toy data set (reads that align to first 1000 bp of *E. coli*). To try SPAdes on this data set, run:

``` bash

    <spades installation dir>/spades.py --test
```

If you added SPAdes installation directory to the `PATH` variable, you can run:

``` bash

    spades.py --test
```

For the simplicity we further assume that SPAdes installation directory is added to the `PATH` variable.

If the installation is successful, you will find the following information at the end of the log:

``` plain

===== Assembling finished. Used k-mer sizes: 21, 33, 55

 * Corrected reads are in spades_test/corrected/
 * Assembled contigs are in spades_test/contigs.fasta
 * Assembled scaffolds are in spades_test/scaffolds.fasta
 * Assembly graph is in spades_test/assembly_graph.fastg
 * Assembly graph in GFA format is in spades_test/assembly_graph.gfa
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
## SPAdes input

SPAdes takes as input paired-end reads, mate-pairs and single (unpaired) reads in FASTA and FASTQ. For IonTorrent data SPAdes also supports unpaired reads in unmapped BAM format (like the one produced by Torrent Server). However, in order to run read error correction, reads should be in FASTQ or BAM format. Sanger, Oxford Nanopore and PacBio CLR reads can be provided in both formats since SPAdes does not run error correction for these types of data.

To run SPAdes 3.13.0 you need at least one library of the following types:

-   Illumina paired-end/high-quality mate-pairs/unpaired reads
-   IonTorrent paired-end/high-quality mate-pairs/unpaired reads
-   PacBio CCS reads

Illumina and IonTorrent libraries should not be assembled together. All other types of input data are compatible. SPAdes should not be used if only PacBio CLR, Oxford Nanopore, Sanger reads or additional contigs are available.

SPAdes supports mate-pair only assembly. However, we recommend to use only high-quality mate-pair libraries in this case (e.g. that do not have a paired-end part). We tested mate-pair only pipeline using Illumina Nextera mate-pairs. See more [here](#hqmp).

Current version SPAdes also supports Lucigen NxSeq® Long Mate Pair libraries, which always have forward-reverse orientation. If you wish to use Lucigen NxSeq® Long Mate Pair reads, you will need Python [regex library](https://pypi.python.org/pypi/regex) to be pre-installed on your machine. You can install it with Python [pip-installer](http://www.pip-installer.org/):

``` bash

    pip install regex
```

or with the [Easy Install](http://peak.telecommunity.com/DevCenter/EasyInstall) Python module:

``` bash

    easy_install regex
```

Notes:

-   It is strongly suggested to provide multiple paired-end and mate-pair libraries according to their insert size (from smallest to longest).
-   It is not recommended to run SPAdes on PacBio reads with low coverage (less than 5).
-   We suggest not to run SPAdes on PacBio reads for large genomes.
-   SPAdes accepts gzip-compressed files.

### Read-pair libraries

By using command line interface, you can specify up to nine different paired-end libraries, up to nine mate-pair libraries and also up to nine high-quality mate-pair ones. If you wish to use more, you can use [YAML data set file](#yaml). We further refer to paired-end and mate-pair libraries simply as to read-pair libraries.

By default, SPAdes assumes that paired-end and high-quality mate-pair reads have forward-reverse (fr) orientation and usual mate-pairs have reverse-forward (rf) orientation. However, different orientations can be set for any library by using SPAdes options.

To distinguish reads in pairs we refer to them as left and right reads. For forward-reverse orientation, the forward reads correspond to the left reads and the reverse reads, to the right. Similarly, in reverse-forward orientation left and right reads correspond to reverse and forward reads, respectively, etc.

Each read-pair library can be stored in several files or several pairs of files. Paired reads can be organized in two different ways:

-   In file pairs. In this case left and right reads are placed in different files and go in the same order in respective files.
-   In interleaved files. In this case, the reads are interlaced, so that each right read goes after the corresponding paired left read.

For example, Illumina produces paired-end reads in two files: `R1.fastq` and `R2.fastq`. If you choose to store reads in file pairs make sure that for every read from `R1.fastq` the corresponding paired read from `R2.fastq` is placed in the respective paired file on the same line number. If you choose to use interleaved files, every read from `R1.fastq` should be followed by the corresponding paired read from `R2.fastq`.

If adapter and/or quality trimming software has been used prior to assembly, files with the orphan reads can be provided as "single read files" for the corresponding read-pair library.

<a name="merged"></a>
If you have merged some of the reads from your paired-end (not mate-pair or high-quality mate-pair) library (using tools s.a. [BBMerge](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or [STORM](https://bitbucket.org/yaoornl/align_test/overview)), you should provide the file with resulting reads as a "merged read file" for the corresponding library.
Note that non-empty files with the remaining unmerged left/right reads (separate or interlaced) must be provided for the same library (for SPAdes to correctly detect the original read length).

In an unlikely case some of the reads from your mate-pair (or high-quality mate-pair) library are "merged", you should provide the resulting reads as a SEPARATE single-read library.

### Unpaired (single-read) libraries

By using command line interface, you can specify up to nine different single-read libraries. To input more libraries, you can use [YAML data set file](#yaml).

Single librairies are assumed to have high quality and a reasonable coverage. For example, you can provide PacBio CCS reads as a single-read library.

Note, that you should not specify PacBio CLR, Sanger reads or additional contigs as single-read libraries, each of them has a separate [option](#inputdata). []()

<a name="pacbio"></a>
### PacBio and Oxford Nanopore reads

SPAdes can take as an input an unlimited number of PacBio and Oxford Nanopore libraries.

PacBio CLR and Oxford Nanopore reads are used for hybrid assemblies (e.g. with Illumina or IonTorrent). There is no need to pre-correct this kind of data. SPAdes will use PacBio CLR and Oxford Nanopore reads for gap closure and repeat resolution.

For PacBio you just need to have filtered subreads in FASTQ/FASTA format. Provide these filtered subreads using `--pacbio` option. Oxford Nanopore reads are provided with `--nanopore` option.

PacBio CCS/Reads of Insert reads or pre-corrected (using third-party software) PacBio CLR / Oxford Nanopore reads can be simply provided as single reads to SPAdes.

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

Note that we assume that SPAdes installation directory is added to the `PATH` variable (provide full path to SPAdes executable otherwise: `<spades installation dir>/spades.py`). []()

<a name="basicopt"></a>
### Basic options

`-o <output_dir> `
    Specify the output directory. Required option.

[]()

<a name="sc"></a>
`--sc `
    This flag is required for MDA (single-cell) data.

[]()

<a name="meta"></a>
`--meta `   (same as `metaspades.py`)
    This flag is recommended when assembling metagenomic data sets (runs metaSPAdes, see [paper](https://genome.cshlp.org/content/27/5/824.short) for more details). Currently metaSPAdes supports only a **_single_** short-read library which has to be **_paired-end_** (we hope to remove this restriction soon). In addition, you can provide long reads (e.g. using `--pacbio` or `--nanopore` options), but hybrid assembly for metagenomes remains an experimental pipeline and optimal performance is not guaranteed. It does not support [careful mode](#correctoropt) (mismatch correction is not available). In addition, you cannot specify coverage cutoff for metaSPAdes. Note that metaSPAdes might be very sensitive to presence of the technical sequences remaining in the data (most notably adapter readthroughs), please run quality control and pre-process your data accordingly.

[]()

<a name="plasmid"></a>
`--plasmid `   (same as `plasmidspades.py`)
    This flag is required when assembling only plasmids from WGS data sets (runs plasmidSPAdes, see [paper](http://biorxiv.org/content/early/2016/04/20/048942) for the algorithm details). Note, that plasmidSPAdes is not compatible with [metaSPAdes](#meta) and [single-cell mode](#sc). Additionally, we do not recommend to run plasmidSPAdes on more than one library. See [section 3.6](#sec3.6) for plasmidSPAdes output details.

[]()

<a name="rna"></a>
`--rna `   (same as `rnaspades.py`)
    This flag should be used when assembling RNA-Seq data sets (runs rnaSPAdes). To learn more, see [rnaSPAdes manual](rnaspades_manual.html).

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
    Tries to reduce the number of mismatches and short indels. Also runs MismatchCorrector – a post processing tool, which uses [BWA](http://bio-bwa.sourceforge.net) tool (comes with SPAdes). This option is recommended only for assembly of small genomes. We strongly recommend not to use it for large and medium-size eukaryotic genomes. Note, that this options is is not supported by metaSPAdes and rnaSPAdes. 

`--continue`
    Continues SPAdes run from the specified output folder starting from the last available check-point. Check-points are made after:

-   error correction module is finished
-   iteration for each specified K value of assembly module is finished
-   mismatch correction is finished for contigs or scaffolds

For example, if specified K values are 21, 33 and 55 and SPAdes was stopped or crashed during assembly stage with K = 55, you can run SPAdes with the `--continue` option specifying the same output directory. SPAdes will continue the run starting from the assembly stage with K = 55. Error correction module and iterations for K equal to 21 and 33 will not be run again. If `--continue` is set, the only allowed option is `-o <output_dir> `.

`--restart-from <check_point>`
    Restart SPAdes run from the specified output folder starting from the specified check-point. Check-points are:

-   `ec` – start from error correction
-   `as` – restart assembly module from the first iteration
-   `k<int>` – restart from the iteration with specified k values, e.g. k55 (not available in RNA-Seq mode)
-   `mc` – restart mismatch correction
-   `last` – restart from the last available check-point (similar to `--continue`)

In contrast to the `--continue` option, you can change some of the options when using `--restart-from`. You can change any option except: all basic options, all options for specifying input data (including `--dataset`), `--only-error-correction` option and `--only-assembler` option. For example, if you ran assembler with k values 21,33,55 without mismatch correction, you can add one more iteration with k=77 and run mismatch correction step by running SPAdes with following options:
`--restart-from k55 -k 21,33,55,77 --mismatch-correction -o <previous_output_dir>`.
Since all files will be overwritten, do not forget to copy your assembly from the previous run if you need it.

`--disable-gzip-output`
    Forces read error correction module not to compress the corrected reads. If this options is not set, corrected reads will be in `*.fastq.gz` format.

[]()

<a name="inputdata"></a>
### Input data

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
    Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) must be provided for the same library for SPAdes to correctly detect the original read length.

`-s <file_name> `
    File with unpaired reads.

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
    If the properties of the library permit, paired reads can be merged using special software.     Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) must be provided for the same library for SPAdes to correctly detect the original read length.

`--pe<#>-s <file_name> `
    File with unpaired reads from paired-end library number `<#>` (`<#>` = 1,2,..,9)
    For example, paired reads can become unpaired during the error correction procedure.

`--pe<#>-<or> `
    Orientation of reads for paired-end library number `<#>` (`<#>` = 1,2,..,9; `<or>` = "fr","rf","ff").
    The default orientation for paired-end libraries is forward-reverse (`--> <--`). For example, to specify reverse-forward orientation for the second paired-end library, you should use the flag: `--pe2-rf `
    Should not be confused with FR and RF strand-specificity for RNA-Seq data (see <a href="rnaspades_manual.html#sec2.3" target="_blank">rnaSPAdes manual</a>). 

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

<a name="lxmp"></a>
**_Lucigen NxSeq® Long Mate Pair libraries_** (see [section 3.1](#sec3.1) for details)

`--nxmate<#>-1 <file_name> `
    File with left reads for Lucigen NxSeq® Long Mate Pair library number `<#>` (`<#>` = 1,2,..,9).

`--nxmate<#>-2 <file_name> `
    File with right reads for Lucigen NxSeq® Long Mate Pair library number `<#>` (`<#>` = 1,2,..,9).

#### Specifying data for hybrid assembly

`--pacbio <file_name> `
    File with PacBio CLR reads. For PacBio CCS reads use `-s` option. More information on PacBio reads is provided in [section 3.1](#pacbio).

`--nanopore <file_name> `
    File with Oxford Nanopore reads.

`--sanger <file_name> `
    File with Sanger reads

`--trusted-contigs <file_name> `
    Reliable contigs of the same genome, which are likely to have no misassemblies and small rate of other errors (e.g. mismatches and indels). This option is not intended for contigs of the related species.

`--untrusted-contigs <file_name> `
    Contigs of the same genome, quality of which is average or unknown. Contigs of poor quality can be used but may introduce errors in the assembly. This option is also not intended for contigs of the related species.

<a name="yaml"></a>
####  Specifying input data with YAML data set file (advanced)

An alternative way to specify an input data set for SPAdes is to create a [YAML](http://www.yaml.org/) data set file. By using a YAML file you can provide an unlimited number of paired-end, mate-pair and unpaired libraries. Basically, YAML data set file is a text file, in which input libraries are provided as a comma-separated list in square brackets. Each library is provided in braces as a comma-separated list of attributes. The following attributes are available:

-   orientation ("fr", "rf", "ff")
-   type ("paired-end", "mate-pairs", "hq-mate-pairs", "single", "pacbio", "nanopore", "sanger", "trusted-contigs", "untrusted-contigs")
-   interlaced reads (comma-separated list of files with interlaced reads)
-   left reads (comma-separated list of files with left reads)
-   right reads (comma-separated list of files with right reads)
-   single reads (comma-separated list of files with single reads or unpaired reads from paired library)
-   merged reads (comma-separated list of files with [merged reads](#merged))

To properly specify a library you should provide its type and at least one file with reads. Orientation is an optional attribute. Its default value is "fr" (forward-reverse) for paired-end libraries and "rf" (reverse-forward) for mate-pair libraries.

The value for each attribute is given after a colon. Comma-separated lists of files should be given in square brackets. For each file you should provide its full path in double quotes. Make sure that files with right reads are given in the same order as corresponding files with left reads.

For example, if you have one paired-end library splitted into two pairs of files:

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
    Read coverage cutoff value. Must be a positive float value, or 'auto', or 'off'. Default value is 'off'. When set to 'auto' SPAdes automatically computes coverage threshold using conservative strategy. Note, that this option is not supported by metaSPAdes.

`--phred-offset <33 or 64>`
    PHRED quality offset for the input reads, can be either 33 or 64. It will be auto-detected if it is not specified.


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

If a single-read library is splitted into several files:

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

The selection of k-mer length is non-trivial for IonTorrent. If the dataset is more or less conventional (good coverage, not high GC, etc), then use our [recommendation for long reads](#sec3.4) (e.g. assemble using k-mer lengths 21,33,55,77,99,127). However, due to increased error rate some changes of k-mer lengths (e.g. selection of shorter ones) may be required. For example, if you ran SPAdes with k-mer lengths 21,33,55,77 and then decided to assemble the same data set using more iterations and larger values of K, you can run SPAdes once again specifying the same output folder and the following options: `--restart-from k77 -k 21,33,55,77,99,127 --mismatch-correction -o <previous_output_dir>`. Do not forget to copy contigs and scaffolds from the previous run. We're planning to tackle issue of selecting k-mer lengths for IonTorrent reads in next versions.

You may need no error correction for Hi-Q enzyme at all. However, we suggest trying to assemble your data with and without error correction and select the best variant.

For non-trivial datasets (e.g. with high GC, low or uneven coverage) we suggest to enable single-cell mode (setting `--sc` option) and use k-mer lengths of 21,33,55. []()

<a name="sec3.4"></a>
## Assembling long Illumina paired reads (2x150 and 2x250)

Recent advances in DNA sequencing technology have led to a rapid increase in read length. Nowadays, it is a common situation to have a data set consisting of 2x150 or 2x250 paired-end reads produced by Illumina MiSeq or HiSeq2500. However, the use of longer reads alone will not automatically improve assembly quality. An assembler that can properly take advantage of them is needed.

SPAdes' use of iterative k-mer lengths allows benefiting from the full potential of the long paired-end reads. Currently one has to set the assembler options up manually, but we plan to incorporate automatic calculation of necessary options soon.

Please note that in addition to the read length, the insert length also matters a lot. It is not recommended to sequence a 300bp fragment with a pair of 250bp reads. We suggest using 350-500 bp fragments with 2x150 reads and 550-700 bp fragments with 2x250 reads.

### Multi-cell data set with read length 2x150

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

### Multi-cell data set with read lengths 2 x 250

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

### Single-cell data set with read lengths 2 x 150 or 2 x 250

The default k-mer lengths are recommended. For single-cell data sets SPAdes selects k-mer sizes 21, 33 and 55.

However, it might be tricky to fully utilize the advantages of long reads you have. Consider contacting us for more information and to discuss assembly strategy.
[]()

<a name="sec3.5"></a>
## SPAdes output

SPAdes stores all output files in `<output_dir> `, which is set by the user.

-   `<output_dir>/corrected/` directory contains reads corrected by BayesHammer in `*.fastq.gz` files; if compression is disabled, reads are stored in uncompressed `*.fastq` files
-   `<output_dir>/scaffolds.fasta` contains resulting scaffolds (recommended for use as resulting sequences)
-   `<output_dir>/contigs.fasta` contains resulting contigs
-   `<output_dir>/assembly_graph.gfa` contains SPAdes assembly graph and scaffolds paths in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md)
-   `<output_dir>/assembly_graph.fastg` contains SPAdes assembly graph in [FASTG format](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf)
-   `<output_dir>/contigs.paths` contains paths in the assembly graph corresponding to contigs.fasta (see details below)
-   `<output_dir>/scaffolds.paths` contains paths in the assembly graph corresponding to scaffolds.fasta (see details below)

Contigs/scaffolds names in SPAdes output FASTA files have the following format:
`>NODE_3_length_237403_cov_243.207`
Here `3` is the number of the contig/scaffold, `237403` is the sequence length in nucleotides and `243.207` is the k-mer coverage for the last (largest) k value used. Note that the k-mer coverage is always lower than the read (per-base) coverage.

In general, SPAdes uses two techniques for joining contigs into scaffolds. First one relies on read pairs and tries to estimate the size of the gap separating contigs. The second one relies on the assembly graph: e.g. if two contigs are separated by a complex tandem repeat, that cannot be resolved exactly, contigs are joined into scaffold with a fixed gap size of 100 bp. Contigs produced by SPAdes do not contain N symbols.

To view FASTG and GFA files we recommend to use [Bandage visualization tool](http://rrwick.github.io/Bandage/). Note that sequences stored in `assembly_graph.fastg` correspond to contigs before repeat resolution (edges of the assembly graph). Paths corresponding to contigs after repeat resolution (scaffolding) are stored in `contigs.paths` (`scaffolds.paths`) in the format accepted by Bandage (see [Bandage wiki](https://github.com/rrwick/Bandage/wiki/Graph-paths) for details). The example is given below.

Let the contig with the name `NODE_5_length_100000_cov_215.651` consist of the following edges of the assembly graph:

``` plain
    >EDGE_2_length_33280_cov_199.702
    >EDGE_5_length_84_cov_321.414'
    >EDGE_3_length_111_cov_175.304
    >EDGE_5_length_84_cov_321.414'
    >EDGE_4_length_66661_cov_223.548
```

Then, `contigs.paths` will contain the following record:

``` plain
    NODE_5_length_100000_cov_215.651
    2+,5-,3+,5-,4+
```


Since the current version of Bandage does not accept paths with gaps, paths corresponding contigs/scaffolds jumping over a gap in the assembly graph are splitted by semicolon at the gap positions. For example, the following record

``` plain
    NODE_3_length_237403_cov_243.207
    21-,17-,15+,17-,16+;
    31+,23-,22+,23-,4-
```

states that `NODE_3_length_237403_cov_243.207` corresponds to the path with 10 edges, but jumps over a gap between edges `EDGE_16_length_21503_cov_482.709` and `EDGE_31_length_140767_cov_220.239`.

The full list of `<output_dir>` content is presented below:

- scaffolds.fasta – resulting scaffolds (recommended for use as resulting sequences)
- contigs.fasta – resulting contigs
- assembly_graph.fastg – assembly graph
- contigs.paths – contigs paths in the assembly graph
- scaffolds.paths – scaffolds paths in the assembly graph
- before_rr.fasta – contigs before repeat resolution

- corrected/ – files from read error correction
    - configs/ – configuration files for read error correction
    - corrected.yaml – internal configuration file
    - Output files with corrected reads

- params.txt – information about SPAdes parameters in this run
- spades.log – SPAdes log
- dataset.info – internal configuration file
- input_dataset.yaml – internal YAML data set file
- K<##>/ – directory containing intermediate files from the run with K=<##>. These files should not be used as assembly results; use resulting contigs/scaffolds in files mentioned above.


SPAdes will overwrite these files and directories if they exist in the specified `<output_dir>`. []()

<a name="sec3.6"></a>
## plasmidSPAdes output

plasmidSPAdes outputs only DNA sequences from putative plasmids. Output file names and formats remain the same as in SPAdes (see [previous](#sec3.5) section), with the following difference. For all contig names in `contigs.fasta`, `scaffolds.fasta` and `assembly_graph.fastg` we append suffix `_component_X`, where `X` is the id of the putative plasmid, which the contig belongs to. Note that plasmidSPAdes may not be able to separate similar plasmids and thus their contigs may appear with the same id. []()

<a name="sec3.7"></a>
## Assembly evaluation

[QUAST](http://cab.spbu.ru/software/quast/) may be used to generate summary statistics (N50, maximum contig length, GC %, \# genes found in a reference list or with built-in gene finding tools, etc.) for a single assembly. It may also be used to compare statistics for multiple assemblies of the same data set (e.g., SPAdes run with different parameters, or several different assemblers).
[]()


<a name="sec4"></a>
# Stand-alone binaries released within SPAdes package

<a name="sec4.1"></a>
## k-mer counting

To provide input data to SPAdes k-mer counting tool `spades-kmercounter ` you may just specify files in [SPAdes-supported formats](#sec3.1) without any flags (after all options) or provide dataset description file in [YAML format](#yaml).

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
## Graph construction
Graph construction tool `spades-gbuilder ` has two mandatory options: dataset description file in [YAML format](#yaml) and an output file name.

Synopsis: `spades-gbuilder <dataset description (in YAML)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>] [-b <value>] [-unitigs|-fastg|-gfa|-spades]`

Additional options are:

`-k <int> `
    k-mer length used for construction (must be odd)

`-t <int> `
    number of threads

`-tmpdir <dir_name>  `
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


<a name="sec4.3"></a>
## Long read to graph aligner
A tool for aligning long reads to the graph `spades-gmapper ` has three mandatory options: dataset description file in [YAML format](#yaml), graph file in GFA format and an output file name.

Synopsis: `spades-gmapper <dataset description (in YAML)> <graph (in GFA)> <output filename> [-k <value>] [-t <value>] [-tmpdir <dir>]`

Additional options are:

`-k <int> `
    k-mer length that was used for graph construction

`-t <int> `
    number of threads

`-tmpdir <dir_name>  `
    scratch directory to use


<a name="sec5"></a>
# Citation
If you use SPAdes in your research, please include [Nurk, Bankevich et al., 2013](http://link.springer.com/chapter/10.1007%2F978-3-642-37195-0_13) in your reference list. You may also add [Bankevich, Nurk et al., 2012](http://online.liebertpub.com/doi/abs/10.1089/cmb.2012.0021) instead.

In case you perform hybrid assembly ussing  PacBio or Nanopore reads, you may also cite [Antipov et al., 2015](http://bioinformatics.oxfordjournals.org/content/early/2015/11/20/bioinformatics.btv688.short). 

If you use multiple paired-end and/or mate-pair libraries you may also cite papers describing SPAdes repeat resolution algorithms [Prjibelski et al., 2014](http://bioinformatics.oxfordjournals.org/content/30/12/i293.short) and [Vasilinetc et al., 2015](http://bioinformatics.oxfordjournals.org/content/31/20/3262.abstract). 

If you use plasmidSPAdes please cite [Antipov et al., 2016](http://biorxiv.org/content/early/2016/04/20/048942).

For rnaSPAdes citation use [Bushmanova et al., 2018](https://www.biorxiv.org/content/early/2018/09/18/420208)</a>.

In addition, we would like to list your publications that use our software on our website. Please email the reference, the name of your lab, department and institution to <spades.support@cab.spbu.ru>.
[]()

<a name="sec6"></a>
# Feedback and bug reports

Your comments, bug reports, and suggestions are very welcomed. They will help us to further improve SPAdes. If you have any troubles running SPAdes, please send us `params.txt` and `spades.log` from the directory `<output_dir>`.

You can leave your comments and bug reports at [our GitHub repository tracker](https://github.com/ablab/spades/issues) or sent it via e-mail: <spades.support@cab.spbu.ru>.

