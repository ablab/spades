<a name="sec1"></a>
# About

This is a support software release for cloudSPAdes publication.

cloudSPAdes is a module of the SPAdes assembler (<http://cab.spbu.ru/software/spades/>) aimed at the data generated using Synthetic Long Read technologies. Currently cloudSPAdes support 10X Genomics Chromium and TELL-seq technologies. This document describes the part of the SPAdes interface which is relevant to the cloudSPAdes module. Please address the main SPAdes manual <http://cab.spbu.ru/files/release3.13.0/manual.html> for details. []()

<a name="sec2"></a>
# Installation

<a name="sec2.1"></a>
# Compiling SPAdes source code

You will need the following libraries to be pre-installed:

-   g++ (version 5.3.1 or higher)
-   cmake (version 2.8.12 or higher)
-   zlib
-   libbz2

If you meet these requirements, you can unpack and build the cloudSPAdes source code:

``` bash

    tar -xzf cloudSPAdes-ismb.tar.gz
    cd cloudSPAdes-ismb
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

which will install SPAdes into `/usr/local/bin`. []()

<a name="sec2.2"></a>
# Verifying your installation 
For testing purposes, SPAdes comes with a toy data set (reads that align to first 1000 bp of *E. coli*). To try SPAdes on this data set, run:

``` bash

    <spades installation dir>/spades.py --test
```

If you added SPAdes installation directory to the `PATH` variable, you can run:

``` bash

    spades.py --test
```

To try SPAdes on a toy SLR dataset, run:

``` bash

    cd <spades installation dir>
    ./spades.py --gemcode1-12 test_dataset_cloudspades/test.fastq  -o cloudspades_test
```
[]()

<a name="sec3"></a>
# Running cloudSPAdes 

To pass SLR reads to this version of SPAdes, please use spades.py with --gemcode1-1 and --gemcode1-2 options for left and right reads respectively (instead of --pe1-1 and --pe1-2 for standard Illumina paired-end reads). Reads should be in FASTQ format with barcodes attached as BC:Z or BX:Z tags:

``` bash

@COOPER:77:HCYNTBBXX:1:1216:22343:0 BX:Z:AAAAAAAAAACATAGT
CCAGGTAGGATTATGGAATTGGTATAAGCGATCAAACTCAATATTTTTGGTGCGGTGACAGACGCCTTCTGGCAGATGATGGGCTTGTCGTAAGTGTGGT
+
GGAGGGAAGGGGIGIIAGAGAGGGGGIAGGGGGGGAGGGGGGGGGGGGAAAGGAGGGGGIGIGGGGGGGAGGAGGIGAIAGGIGGGGIGGGGGGGGGGGG
```

Here is the full list of options related to SLR library processing:

`--gemcode<#>-12 <file_name> `
    File with interlaced reads for SLR library number `<#>` (`<#>` = 1,2,..,9). 

`--gemcode<#>-1 <file_name> `
    File with left reads for SLR library number `<#>` (`<#>` = 1,2,..,9).

`--gemcode<#>-2 <file_name> `
    File with right reads for SLR library number `<#>` (`<#>` = 1,2,..,9).

Otherwise, cloudSPAdes has the same command line interface as SPAdes.

[]()