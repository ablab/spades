<!-- dx-header -->
# SPAdes Assembler (DNAnexus Platform App)

Summary

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
http://wiki.dnanexus.com/.
<!-- /dx-header -->

SPAdes - St. Petersburg Genome Assembler - is intended for both standard
isolates and single-cell MDA bacteria assemblies.

## Supported data types

Current version of SPAdes works only with Illumina reads. Support for other
technologies (e.g. Roche 454, IonTorrent, PacBio) is currently in progress and
probably will be included in one of the next releases.

SPAdes supports paired-end reads as well as unpaired reads. So far SPAdes takes
as input only one paired-end library. We are currently working on mate-pairs and
multiple libraries support – it is likely to come in the next major release.

Also note, that SPAdes was initially designed for single-cell and standard
bacterial data sets and is not intended for larger genomes (e.g. mammalian size
genomes) and metagenomic projects. For such purposes you can use it at your own
risk.

Consult the following application note
http://bioinf.spbau.ru/en/spades-rl250-assembly on assembling long (2x150 and
2x250) Illumina reads.

## SPAdes pipeline

SPAdes comes in three separate modules:

* BayesHammer – read error correction tool, which works well on both single-cell
  and standard data sets.
* SPAdes – iterative short-read genome assembly module; by default consecutively
  iterates through K value equal to 21, 33 and 55.
* MismatchCorrector – a tool which improves mismatch and short indel rates in
  resulting contigs and scaffolds; this module uses BWA tool
  [Li H. and Durbin R., 2009] and is turned on by default.

We recommend to run SPAdes with BayesHammer to obtain high-quality single-cell
assemblies. However, if you use your own read correction tool (e.g. Quake) it is
possible to turn BayesHammer off.

## Inputs

SPAdes takes as input forward-reverse paired-end reads as well as single
(unpaired) reads in FASTA or FASTQ format. However, in order to run read error
correction, reads should be in FASTQ format. Currently SPAdes accepts only one
paired-end library.

The input dataset can be specified using the following options:

* **Left paired end reads** ``left_reads``: ``file`` -  File(s) with forward reads.
* **Right paired end reads** ``right_reads``: ``file`` - File(s) with reverse reads.
* **Unpaired reads** ``single_reads``: ``file`` - File(s) with unpaired reads.
* **Single Cell Dataset** ``is_single_cell``: ``boolean`` - This flag is
  required for MDA (single-cell) data.

## Pipeline options

* **Run assembler only** ``is_only_assembler``: ``boolean`` - Runs assembly
  module only without read error correction. Use for corrected reads.
* **Run assembler in careful mode** ``is_careful_mode``: ``boolean`` - Tries to
  reduce number of mismatches and short indels. Also runs MismatchCorrector – a
  post processing tool, which uses BWA tool (comes with SPAdes).

## Outputs

* **Contigs** ``contigs``: ``file`` - Resulting contigs.
* **Scaffolds** ``scaffolds``: ``file`` - Resulting scaffolds.
