# Quick start

- SPAdes is an assembler for second-generation sequencing data (Illumina or IonTorrent). PacBio and Nanopore reads are supported *only* as supplementary data. SPAdes can assemble genomes, metagenomes, transcriptomes, viral genomes etc.

- Download SPAdes binaries for [Linux](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz) or [MacOS](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Darwin.tar.gz). You can also compile SPAdes from [source](https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5.tar.gz) (requires g++ 9.0+, cmake 3.16+, zlib and libbz2). SPAdes requires only Python 3.8+ to be installed.

- Test your SPAdes installation by running

```
    bin/spades.py --test
```

- A single paired-end library (separate files, gzipped):

```
    bin/spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- A single paired-end library (interlaced reads):

```
    bin/spades.py --12 interlaced.fastq -o output_folder
```

- Two paired-end libraries (separate files):

```
    bin/spades.py --pe1-1 1_left.fastq --pe1-2 1_right.fastq --pe2-1 2_left.fastq --pe2-2 2_right.fastq -o output_folder
```

- IonTorrent data:
```
    bin/spades.py --iontorrent -s it_reads.fastq -o output_folder
```

- A paired-end library coupled with long PacBio reads:

```
    bin/spades.py -1 left.fastq.gz -2 right.fastq.gz --pacbio pb.fastq -o output_folder
```

- Available assembly modes: `--isolate`, `--sc`, `--plasmid`, `--meta`, `--metaplasmid`, `--metaviral`, `--rna`, `--rnaviral`, `--corona`, `--bio`.




