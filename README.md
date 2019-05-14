# cloudSPAdes: assembly of synthetic long reads using de Bruijn graphs

## MANUAL

### Overview

cloudSPAdes is a module of the SPAdes assembler (<http://cab.spbu.ru/software/spades/>) aimed at genome assembly from the data generated using Synthetic Long Read (SLR) technologies. This software release consists of a main SPAdes executable script for genome assembly (assembler/spades.py) with additional command line options allowing to assemble SLR datasets. 

### Input

The tool supports SLR libraries produced using 10X Genomics Chromium and TELL-seq technologies. To pass an SLR library to cloudSPAdes, please use spades.py with `--gemcode1-1` and `--gemcode1-2` options for left and right reads respectively (instead of `--pe1-1` and `--pe1-2` for standard Illumina paired-end reads). SLR library should be in FASTQ format with barcodes attached as BC:Z or BX:Z tags:

``` 
@COOPER:77:HCYNTBBXX:1:1216:22343:0 BX:Z:AAAAAAAAAACATAGT
CCAGGTAGGATTATGGAATTGGTATAAGCGATCAAACTCAATATTTTTGGTGCGGTGACAGACGCCTTCTGGCAGATGATGGGCTTGTCGTAAGTGTGGT
+
GGAGGGAAGGGGIGIIAGAGAGGGGGIAGGGGGGGAGGGGGGGGGGGGAAAGGAGGGGGIGIGGGGGGGAGGAGGIGAIAGGIGGGGIGGGGGGGGGGGG
```

Currently cloudSPAdes supports only a single SLR library. Please use regular SPAdes for any input libraries combination that does not include SLR libraries.

#### Command line options

Here is the full list of options related to SLR library processing:

`--gemcode<#>-12 <file_name> `
    File with interlaced reads for SLR library number `<#>` (`<#>` = 1,2,..,9). 

`--gemcode<#>-1 <file_name> `
    File with left reads for SLR library number `<#>` (`<#>` = 1,2,..,9).

`--gemcode<#>-2 <file_name> `
    File with right reads for SLR library number `<#>` (`<#>` = 1,2,..,9).

Otherwise, cloudSPAdes has the same command line interface as SPAdes.

### Output

cloudSPAdes produces same output files as the regular SPAdes. Please address <http://cab.spbu.ru/files/release3.13.1/manual.html> for details.

### Installation

Installation from the source code is the same as in the regular SPAdes. 

### Example

To try cloudSPAdes on a toy SLR dataset, run:

``` 
    cd <cloudSPAdes installation dir>
    ./spades.py --gemcode1-12 test_dataset_cloudspades/test.fastq  -o cloudspades_test
```
