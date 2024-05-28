# SPAdes Genome Assembler changelog

## SPAdes 4.0.0, 3 June 2024
- Python3.8 is now a minimal required version, python2 is deprecated;
- SPAdes now supports SRA files as input;
- PathRacer is now released as a part of SPAdes package;
- [BinSPreader](https://ablab.github.io/spades/binspreader.html) is also along with SPAdes package;
- New gfa-split tool;
- New [wasterwaterSPAdes](https://ablab.github.io/spades/running.html#-sewage) mode;
- Native support for running on Apple Silicon processors;
- Native support for running  on aarch64/linux;
- Binaries are now built using ManyLinux 2.28;
- Parallel & improved GapCloser algorithm;
- Faster graph simplification;
- Allow external projects integration with SPAdes;
- Use SPOA for long reads consensus;
- SPAdes outputs assembly graph in GFA v1.2 format by default, older v1.1 format is also available;
- SPAdes tags circular paths in assembly graph with `TP:Z:circular` flag;
- SPAdes now uses `zlib-ng` for faster gzip decompression of input files;
- Fixed compilation with gcc 13;
- Removed deprecated truSPAdes functionality;
- Reworked and improved documentation;
- Lots of other fixes and improvements.


## SPAdes 3.15.5, 14 July 2022

- NEW: Support DNA HMMs.

- FIX: use-after-free in PB aligner.

- FIX: duplicated sequences in metaplasmid / metaviral mode.


## SPAdes 3.15.4, 1 February 2022

- FIX: MacOS Monterey memory limit failure.

- FIX: upgrade pyyaml to run correctly with Python 3.10.

- FIX: WSL for py2 check.

- FIX: A few stability fixes.


## SPAdes 3.15.3, 22 July 2021

- FIX: trusted contigs failure.

- FIX: clarification & refining the output of bgcSPAdes and coronaSPAdes.

- FIX: usage of >9 libraries in a single SPAdes run.

- FIX: improvements in `spades-read-filter` tool.


## SPAdes 3.15.2, 8 March 2021

- FIX: meta-viral pipeline bugs.

- FIX: coronaspades.py wrapper, copy proper files to the output folder.

- FIX: coronaSPAdes instability


## SPAdes 3.15.1, 18 February 2021

- FIX: Gap closer failure when using multiple libraries.

- FIX: Gap closer excessive memory consumption.

- IMPROVE: coronaSPAdes output.


## SPAdes 3.15.0, 11 January 2021

- NEW: CoronaSPAdes pipeline for assembly of full-length coronaviridae genomes from the transcriptomic and metatranscriptomic data.

- NEW: Meta-Viral and RNA-Viral pipelines for identifying viral genomes for metagenomic and metatranscriptomic data.

- NEW: Novel algorithm for trusted contig usage.

- NEW: Switched to [mimalloc](https://github.com/microsoft/mimalloc) memory allocator.

- NEW: PlasmidSPAdes and bgcSPAdes now support assembly graph as an input.

- CHANGE: Significant improvements and fixes for metaplasmid pipeline.

- CHANGE: Multiple performance improvements in simplification and repeat resolving procedures.

- DEPRECATED: Support for Lucigen NxSeq® Long Mate Pair reads.

- DEPRECATED: truSPAdes pipeline for TruSeq barcode assembly (still present in this release but no longer supported).


## SPAdes 3.14.1, 27 April 2020

- FIX: metaplasmidSPAdes contig output.

- FIX: read filtering binary.

- FIX: biosyntheticSPAdes pipeline fixed.

- FIX: fixed truSPAdes for Python 3.6+.

- FIX: bug in the internal mismatch correction procedure.

- FIX: Soft and hard-filtered transcripts are now copied to the output folder in rnaSPAdes.

- FIX: Several usability fixes in `spades.py'.

- FIX: meta-plasmid options added to the manual.

- FIX: several minor fixes in the user manual.


## SPAdes 3.14.0, 27 December 2019

- NEW: BiosyntheticSPAdes pipeline for predicting Biosynthetic Gene Clusters.

- NEW: Hybrid transcriptome assembly with rnaSPAdes.

- NEW: Plasmid detection from metagenomic samples.

- NEW: Special `--isolate` option for assembly of standard datasets with good coverage (>100x).

- NEW: Standalone tool for reads filtration based on k-mer coverage.

- NEW: Standalone tool for estimating approximate number of unique k-mers in reads.

- CHANGE: Improved SPAligner tool.

- CHANGE: Reworked python code, faster sequence transfer between different k-mer stages.

- CHANGE: Multiple performance improvements in graph construction and simplification procedures.

- FIX: BWA aligner failure for large graphs.

- FIX: Failure when additional paired-end libraries with reads shorter than final k-mer length are provided.


## SPAdes 3.13.2, 31 October 2019

- FIX: Incorrect k-mer size estimation in rnaSPAdes.


## SPAdes 3.13.1, 11 April 2019

- CHANGE: Removed BayesHammer from rnaSPAdes pipeline.

- CHANGE: Improved rnaSPAdes performance on large datasets.

- FIX: Failure in contig output in rnaSPAdes.


## SPAdes 3.13.0, 11 October 2018

- CHANGE: Switched to multi-k-mer mode in rnaSPAdes, k-mer values are detected automatically based on read length.

- CHANGE: Added manual as README.md in markdown format for github.

- FIX: Updated BWA and switched to RopeBWT, which allows to handle large graphs.

- FIX: Assert for path.length() > 0 in rnaSPAdes.

- FIX: CQF and k-mer counting.


## SPAdes 3.12.0, 14 May 2018

- NEW: Support for merged paired-end reads.

- NEW: Experimental pipeline for metagenome hybrid assemblies.

- NEW: Standalone graph builder application.

- NEW: Standalone k-mer counting application.

- NEW: Standalone long read to graph aligner.

- CHANGE: Significant improvements in hybrid assembly pipeline.

- CHANGE: Faster read alignment using BWA.

- CHANGE: Improvements in metaSPAdes results.

- CHANGE: More sensitive results for rnaSPAdes.

- CHANGE: All binaries for SPAdes pipeline steps now have `spades-` prefix in its name.

- CHANGE: Better running time and RAM consumption for graph construction stage.

- CHANGE: Overall performance improvements.

- FIX: K value estimation for rnaSPAdes.

- DEPRECATED: dipSPAdes pipeline for highly polymorphic diploid genomes (still present in the release but no longer supported).


## SPAdes 3.11.1, 1 October 2017

- FIX: Handling spaces in path during mismatch correction.

- FIX: Python3 support in rnaSPAdes.

- FIX: K value estimation for long reads.

- FIX: Processing long reads alignments.


## SPAdes 3.11.0, 1 September 2017

- NEW: Support for strand-specific RNA-Seq data in rnaSPAdes.

- NEW: Coverage based isoform detection in rnaSPAdes.

- NEW: Reworked IonHammer read error correction module.

- CHANGE: Improved tandem repeat resolution accuracy.

- CHANGE: Better performance of exSPAnder module.

- CHANGE: metaSPAdes pipeline improvements.

- CHANGE: Better running time and RAM consumption for the entire pipeline.

- FIX: Incomplete paths in GFA output.

- FIX: Mismatch and indel rate in careful mode for isolate datasets (esp. low covered ones).

- FIX: Occasional hanging of edge disconnection procedure in metaSPAdes.


## SPAdes 3.10.1, 1 March 2017

- FIX: Build for MacOS.

- FIX: Minor bugs in hybridSPAdes pipeline.

- FIX: `--continue` option for metaSPAdes.

- FIX: `--tmp-dir` is now works correctly for MismatchCorrector.

- FIX: `Assertion `overlap <= k_' failed` in rnaSPAdes and metaSPAdes.

- FIX: `Assertion `path.Length() > 0' failed` in metaSPAdes.


## SPAdes 3.10.0, 27 January 2017

- NEW: Scaffolding algorithm for mate-pairs and long reads.

- NEW: Contigs and graph output in GFA format.

- CHANGE: Better running time and RAM consumption for all pipelines.

- CHANGE: Improvements in metagenomic pipeline.

- CHANGE: Improved isoform detection algorithm in rnaSPAdes.


## SPAdes 3.9.1, 4 December 2016

- FIX: macOS Sierra crash.

## SPAdes 3.9.0, 23 July 2016

- NEW: rnaSPAdes pipeline for de novo transcriptome assembly from RNA-Seq data.

- CHANGE: Improved memory consumption in metagenomic pipeline.

- FIX: Several minor bugs.


## SPAdes 3.8.2, 10 July 2016

-  FIX: Several minor bug-fixes for metaSPAdes and SPAdes pipelines.


## SPAdes 3.8.1, 8 June 2016

- FIX: plasmidSPAdes now works with PacBio/Nanopore reads.


## SPAdes 3.8.0, 1 June 2016

- NEW: Added plasmidSPAdes &ndash;  a pipeline designed for extracting and assembling plasmids from WGS data sets.

- CHANGE: Significant improvements in metaSPAdes performance.

- CHANGE: Improved running time and RAM consumption.


## SPAdes 3.7.1, 8 March 2016

- FIX: MismatchCorrector fixed for MaxOS.


## SPAdes 3.7.0, 24 February 2016

- NEW: metaSPAdes metagenomic pipeline.

- CHANGE: improved performance for both error correction and assembly stages.

- FIX: Multiple bug fixes.


## SPAdes 3.6.2, 20 November 2015

- NEW: Contigs/scaffolds paths for assembly_graph.fastg in Bandage-supported format.

- FIX: Multithreaded MismatchCorrector.

- FIX: BayesHammer bug fixes.

- FIX: Python 3.5 support; python 3 support for truSPAdes.

## SPAdes 3.6.1, 4 October 2015

- CHANGE: No misleading FASTG files, only assembly graph is saved in FASTG format.

- FIX: Multiple bugfixes.

## SPAdes 3.6.0, 17 August 2015

- NEW: Added truSPAdes &ndash; an assembler for short reads produced by Illumina TruSeq Long Read technology.

- CHANGE: Better running time, less RAM consumption and improved results for BayesHammer error correction module.

- CHANGE: Improvements and bugfixes in repeat resolution and scaffolding modules.

- CHANGE: Improvements and bugfixes in dipSPAdes.

- CHANGE: MismatchCorrector now uses bwa-mem.

- FIX: Bugfixes in MismatchCorrector.


## SPAdes 3.5.0, 7 December 2014

- NEW: New MismatchCorrector module.

- NEW: Support for Oxford Nanopore long reads.

- NEW: Support for Lucigen NxMate mate-pair libraries.

- NEW: Possibility to specify coverage cutoff: automatic and manual.

- CHANGE: Better running time.

- CHANGE: Improved RAM consumption.

- CHANGE: High-quality mate-pairs are now assumed to have forward-revers orientation (same as paired-end).

- FIX: Fixed FASTG format.


## SPAdes 3.1.1, 27 May 2014

- FIX: Several improvements in IonHammer.

- FIX: Fixed a few minor bugs in repeat resolution and scaffolding.

## SPAdes 3.1.0, 27 May 2014

- NEW: Mate-pair only assembly with high-quality libraties.

- NEW: Support for BAM files.

- CHANGE: Improved IonTorrent pipeline.

- CHANGE: Better quality and higher performance when using mate-pairs.

- FIX: Fixed dipSPAdes bugs and user interface.


## SPAdes 3.0.0, 29 December 2013

- NEW: Module for assemblying diplod highly polymorphic genomes.

- NEW: Support for PacBio reads.

- NEW: Support for IonTorrent reads.

- NEW: Support for Sanger reads and additional contigs.

- NEW: Possibility to restart SPAdes starting from the specified check-point with the `--restart-from` option.

- NEW: Output contigs/scaffolds in FASTA and FASTG.

- CHANGE: Improved algorithm for mate-pair repeat resolution and scaffolding.

- CHANGE: Improved N50 and misassembly rate for single-cell data sets with low genome fraction.

- FIX: User-friendly handling for errors in mismatch corrector.

- REMOVE: Rectangle graph repeat resolution module.


## SPAdes 2.5.1, 10 September 2013

- NEW: Python 3.2 and 3.3 compatibility.

- NEW: Possibility to continue SPAdes run starting from the last check-point with the `--continue` option.

- CHANGE: Decreased memory consumption for error correction module.

- CHANGE: Improved misassembly rate for single-cell data sets with low genome fraction.

- FIX: User-friendly handling for the case when paired reads do not align to the assembly graph.


## SPAdes 2.5.0, 2 July 2013

- NEW: Multiple paired-end and mate-pair libraries.

- NEW: Recipe for assembling Illumina 2x150bp and 2x250bp reads.

- CHANGE: Improved mismatch and indel rate.


## SPAdes 2.4.0, 26 February 2013

- NEW: Mismatch correction post-processing module.

- NEW: Rectangle graph repeat resolution module as an option.

- NEW: Build for Mac OS.

- CHANGE: Improved assembly quality of standard (isolate) data sets.

- CHANGE: Decreased memory consumption in error correction module (14 Gb instead of 24 Gb on E.coli test dataset).

- REMOVE: SAM-file generation.

## SPAdes 2.3.0, 30 October 2012

- NEW: Generate scaffolds alongside with contigs.

- CHANGE: Use N instead of A, C, G, T for the variations in repeats.

- CHANGE: Memory requirements for E.coli test dataset decreased from 35 Gb to 24 Gb of RAM.

- CHANGE: output_dir is a required command line parameter instead of project_name.

- CHANGE: Simplified output directory structure.

- CHANGE: CMake 2.8 is required instead of 2.6.

- REMOVE: No dependency from boost library.

## SPAdes 2.2.1, 20 August 2012

- FIX: Avoid `Verification of expression 'v1 == conjugate(v2)' failed` error.


## SPAdes 2.2.0, 2 August 2012

- NEW: No special binaries for different K values.

- NEW: Great improvements in error correction tool BayesHammer.

- CHANGE: Memory requirements for E.coli test dataset decreased from 85 Gb to 35 Gb of RAM.

- CHANGE: Only 1 iteration of BayesHammer by default.

- NEW: Improved assembly quality.


## SPAdes 2.1.0, 28 May 2012

- NEW: Support multi-threading.

- NEW: Improved algorithms.

- NEW: Command-line interface.

- CHANGE: Quality assessment separated from the core pipeline.

- REMOVE: No support for debian and RPM packages.


## SPAdes 2.0.1, 26 Apr 2012

- FIX: Quality tool fixed.


## SPAdes 2.0.0, 18 Apr 2012

- Initial release.
