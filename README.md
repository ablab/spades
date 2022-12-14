# SpLitteR - Repeat resolution in assembly graph using synthetic long reads

## MANUAL

### Description

SpLitteR is a tool that uses synthetic long reads (SLRs) to improve the contiguity of HiFi assemblies. Given a SLR library and a HiFi assembly graph in the GFA format, SpLitteR resolves repeats in the assembly graph using linked-reads and generates a simplified (more contiguous) assembly graph with corresponding scaffolds.

### Dependencies

-   g++ (version 5.3.1 or higher)
-   cmake (version 3.12 or higher)
-   zlib
-   libbz2

### Installation

``` 
cd spades/assembler/
mkdir build && cd build && cmake ../src
make splitter
```
Now to run SpLitteR move to folder `assembler/` and execute

`build/bin/splitter`

### Input

#### Format

The tool requires 

- Assembly graph file in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), with scaffolds included as path lines.
- SLR library in YAML format. The tool supports SLR libraries produced using 10X Genomics Chromium and TELL-seq technologies. SLR library should be in FASTQ format with barcodes attached as BC:Z or BX:Z tags:

``` 
@COOPER:77:HCYNTBBXX:1:1216:22343:0 BX:Z:AAAAAAAAAACATAGT
CCAGGTAGGATTATGGAATTGGTATAAGCGATCAAACTCAATATTTTTGGTGCGGTGACAGACGCCTTCTGGCAGATGATGGGCTTGTCGTAAGTGTGGT
+
GGAGGGAAGGGGIGIIAGAGAGGGGGIAGGGGGGGAGGGGGGGGGGGGAAAGGAGGGGGIGIGGGGGGGAGGAGGIGAIAGGIGGGGIGGGGGGGGGGGG
```

For example, if you have an SLR library

``` bash

    lib_slr_1.fastq.gz
    lib_slr_2.fastq.gz
```

YAML file should look like this:

``` bash

    [
      {
        orientation: "fr",
        type: "clouds10x",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_slr_2.fastq.gz" 
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_slr_1.fastq.gz"
        ]
      }
    ]
```

#### Command line

Synopsis: `splitter <graph (in binary or GFA)> <SLR library description (in YAML)> <path to output directory> [OPTION...]`

Main options:

- `-t` Number of threads to use (default: 1/2 of available threads)
- `--mapping-k` k-mer length for read mapping (default: 31)
- `-Gmdbg|-Gblunt` Assembly graph type (mDBG or blunted)
- `-Mdiploid|-Mmeta` Repeat resolution mode (diploid or meta)

Barcode index construction:
- `--count-threshold` Minimum number of reads for barcode index
- `--frame-size` Resolution of the barcode index
- `--length-threshold` Minimum scaffold graph edge length (meta mode option)
- `--linkage-distance` Reads are assigned to the same fragment on long edges based on the linkage distance
- `--min-read-threshold` Minimum number of reads for path cluster extraction
- `--relative-score-threshold` Relative score threshold for path cluster extraction

Repeat resolution:
- `--score` Score threshold for link index.
- `--tail-threshold` Barcodes are assigned to the first and last <tail_threshold> nucleotides of the edge.

Developer options:
- `--ref` Reference path for repeat resolution evaluation
- `--statistics` Produce additional read cloud library statistics
- `--bin-load` Load binary-converted reads from tmpdir
- `--debug` Produce lots of debug data
- `--tmp-dir` Scratch directory to use
- `-h, --help ` Print help message

### Output

SpLitteR stores all output files in output directory `<output_dir> `, which is set by the user.

- `<output_dir>/assembly_graph.gfa` input assembly graph in mDBG encoding
- `<output_dir>/resolved_graph.gfa` output assembly graph after repeat resolution
- `<output_dir>/contigs.fasta` output scaffolds

In addition

- `<output_dir>/edge_transform.tsv` map from input graph edges to resolved graph edges
- `<output_dir>/vertex_stats.tsv` Statistics for complex vertices
- `<output_dir>/resolved_graph.fasta` Sequences of the resolved graph edges
