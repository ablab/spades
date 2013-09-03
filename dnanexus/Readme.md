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

SPAdes supports paired-end reads, mate-pairs and unpaired reads. SPAdes can 
take as input several paired-end and mate-pair libraries simultaneously.

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

SPAdes takes as input paired-end reads, mate-pairs and single (unpaired) reads 
in <a href="https://wiki.dnanexus.com/Types/Reads">DNAnexus Reads</a> format. 
Reads in FASTQ formats can be converted into DNAnexus Reads using 
<a href="https://platform.dnanexus.com/app/reads_importer">FASTQ Reads Importer</a> 
App. Paired-end reads and mate-pairs orientation can be specified during the 
conversion. By default, SPAdes assumes that paired-end reads have 
forward-reverse (FR) orientation and mate-pairs have reverse-forward (RF) 
orientation. SPAdes accepts up to five different paired-end libraries and also 
up to five different mate-pair ones.

The input dataset can be specified using the following options:

* **List of paired-end reads** ``paired_reads``: ``array:gtable`` -  array of DNAnexus Reads with paired-end reads.
* **List of mate-pairs** ``mate_pairs``: ``array:gtable`` -  array of DNAnexus Reads with mate-pairs.
* **List of unpaired reads** ``unpaired_reads``: ``array:gtable`` -  array of DNAnexus Reads with unpaired reads.
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
