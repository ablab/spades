# About SPAdes

![SPAdes](spades.png){ align=right }

SPAdes - St. Petersburg genome assembler - a versatile toolkit designed for assembling and analyzing sequencing data from
Illumina and IonTorrent technologies. In addition, most of SPAdes pipelines support a hybrid mode allowing the use of
long reads (PacBio and Oxford Nanopore) as supplementary data.

SPAdes package provides pipelines for DNA assembly of isolates and single-cell bacteria, as well as of
metagenomic and transcriptomic data. Additional modes allow to recover bacterial plasmids and RNA viruses,
perform HMM-guided assembly and more. SPAdes package also includes supplementary tools for efficient
k-mer counting and k-mer-based read filtering, assembly graph construction and simplification,
sequence-to-graph alignment and metagenomic binning refinement.

SPAdes version {{ spades_version() }} was released under GPLv2 on July 14th, 2022 and can be downloaded [here](https://github.com/ablab/spades/releases/latest/).

The latest SPAdes paper describing various pipelines in a protocol format is available [here](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102).
