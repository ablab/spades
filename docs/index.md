# About SPAdes

![SPAdes](spades.png){ align=right }

SPAdes - St. Petersburg genome assembler - is a versatile toolkit designed for assembly and analysis of sequencing data.
SPAdes is primarily developed for second-generation technologies (Illumina and IonTorrent), but most of its pipelines support hybrid mode, i.e. allow using long reads (PacBio and Oxford Nanopore) as a supplementary data.

SPAdes package contains assembly pipelines for isolated and single-cell bacterial, as well as metagenomic and transcriptomic data.
Additional modes allow to discover bacterial plasmids and RNA viruses, as well as perform HMM-guided assembly.
Besides, SPAdes package includes supplementary tools for efficient k-mer counting and filtering, assembly graph construction and simplification, sequence-to-graph alignment and metagenomic binning refinement.

SPAdes version 3.15.5 was released under GPLv2 on July 14th, 2022 and can be downloaded [here](https://github.com/ablab/spades/releases/latest/).

The latest SPAdes paper describing various pipelines in a protocol format is available [here](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102).
