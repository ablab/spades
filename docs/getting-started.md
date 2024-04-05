# Quick start

SPAdes is an assembler for second-generation sequencing data (Illumina or IonTorrent). PacBio and Nanopore reads are supported *only* as supplementary data. SPAdes can assemble genomes, metagenomes, transcriptomes, viral genomes etc.

1. Download SPAdes binaries for
   [Linux](https://github.com/ablab/spades/releases/latest/) or [MacOS](https://github.com/ablab/spades/releases/latest/). 
   You can also compile SPAdes from [source](installation.md#downloading-and-compiling-spades-source-code) (requires g++ 9.0+, cmake 3.16+, zlib and libbz2). SPAdes requires only Python 3.8+ to be installed.

2. Test your SPAdes installation by running
``` bash
    bin/spades.py --test
```

## Useful one-liners
- A single paired-end library (separate files, gzipped):
``` bash
bin/spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- A single paired-end library (interlaced reads):
``` bash
bin/spades.py --12 interlaced.fastq -o output_folder
```

- Two paired-end libraries (separate files):
``` bash
bin/spades.py --pe1-1 1_left.fastq --pe1-2 1_right.fastq \
              --pe2-1 2_left.fastq --pe2-2 2_right.fastq \
              -o output_folder
```

- A paired-end library coupled with long PacBio reads:
``` bash
bin/spades.py -1 left.fastq.gz -2 right.fastq.gz \
              --pacbio pb.fastq \
              -o output_folder
```

- Assemble a uniformly covered isolate bacterial genome :
``` bash
bin/spades.py --isolate -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- Assemble a metagenome:
``` bash
bin/spades.py --meta -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- Assemble a transcriptome:
``` bash
bin/spades.py --rna -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- Assemble an RNA viral genome:
``` bash
bin/spades.py --rnaviral -1 left.fastq.gz -2 right.fastq.gz -o output_folder
```

- Available assembly modes: `--isolate`, `--sc`, `--plasmid`, `--meta`, `--metaplasmid`, `--metaviral`, `--rna`, `--rnaviral`, `--corona`, `--bio`.




