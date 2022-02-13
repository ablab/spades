# BinSPreader - binning refining using assembly graphs

## MANUAL

### Description

BinSPreader is a novel tool that attempts to refine metagenome-assembled genomes (MAGs) obtained from existing tools. BinSPreader exploits the assembly graph topology and other connectivity information, such as paired-end and Hi-C reads, to refine the existing binning, correct binning errors, propagate binning from longer contigs to shorter contigs and infer contigs belonging to multiple bins.

### Dependencies

-   g++ (version 5.3.1 or higher)
-   cmake (version 3.12 or higher)
-   zlib
-   libbz2

### Installation 

``` 
cd spades/assembler/
mkdir build && cd build && cmake ../src
make bin-refine
```
Now to run BinSPreader move to folder `assembler/` and execute 

`build/bin/bin-refine`

### Input

The tool requires initial binning to refine, as well as assembly graph as a source of information for refining. Optionally, BinSPreader can be provided with multiple Hi-C and/or paired-end libraries.

Required positional arguments: 

- Assembly graph file in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), with scaffolds included as path lines. Alternatively, scaffold paths can be provided separately using `--path` option in the `.paths` format accepted by Bandage (see [Bandage wiki](https://github.com/rrwick/Bandage/wiki/Graph-paths) for details). 
- Binning output from an existing tool (in `.tsv` format)

Synopsis: `bin-refine <graph (in GFA)> <binning (in .tsv)> <output directory> [OPTION...]`

Main options:

- `--paths` provide contigs paths from file separately from GFA
- `--dataset` Dataset in YAML format (see #yaml) describing Hi-C and paired-end reads

- `-l` L Library index (0-based, default: 0). Only the library specified by this index will be used.
- `-t` T # of threads to use (default: 1/2 of available threads)
- `-e` E convergence relative tolerance threshold (default: 1e-5)
- `-n` ITERATIONS maximum number of iterations (default: 5000)
- `-m` allow multiple bin assignment (defalut: false)
- `-Smax|-Smle` simple maximum or maximum likelihood binning assignment strategy (default: max likelihood)
- `-Rcorr|-Rprop` Select propagation or correction mode (default: correction)
- `--cami` use CAMI bioboxes binning format
- `--zero-bin` emit zero bin for unbinned sequences
- `--tall-multi` use tall table for multiple binning result
- `--bin-dist` estimate pairwise bin distance (could be slow on large graphs!)
- `-la` LA labels correction regularization parameter for labeled data (default: 0.6)

Sparse propagation options:
- `--sparse-propagation` Gradually reduce regularization parameter from binned to unbinned edges. Recommended for sparse binnings with low assembly fraction.
- `--no-unbinned-bin` Do not create a special bin for unbinned contigs. More agressive strategy.
- `-ma, --metaalpha` Regularization parameter for sparse propagation procedure. Increase/decrease for more agressive/conservative refining (default: 0.6)
- `-lt, --length-threshold` LENGTH_THRESHOLD Binning will not be propagated to edges longer than threshold
- `-db' --distance-bound` DISTANCE_BOUND Binning will not be propagated further than bound from initially binned edges

Read splitting options:
- `-r, --reads` Split reads according to binning. Can be used for reassembly.
- `-b, --bin-weight` BIN_WEIGHT Reads bin weight threshold (default: 0.1).

Developer options:
- `--bin-load` Load binary-converted reads from tmpdir
- `--debug` produce lots of debug data
- `--tmp-dir` TMP_DIR scratch directory to use
- `-h, --help ` print help message

### BinSPreader output

BinSPreader stores all output files in output directory `<output_dir> `, which is set by the user.

- `<output_dir>/binning.tsv` contains refined binning in `.tsv` format
- `<output_dir>/bin_stats.tsv` contains various per-bin statistics
- `<output_dir>/bin_weights.tsv` contains soft bin weights per contig
- `<output_dir>/edge_weights.tsv` contains soft bin weights per edge

In addition

- `<output_dir>/bin_dist.tsv` contains refined bin distance matrix (if `--bin-dist` was used)
- `<output_dir>/bin_label_1.fastq, <output_dir>/bin_label_2.fastq` read set for bin labeled by `bin_label` (if `--reads` was used)
- `<output_dir>/pe_links.tsv` list of paired-end links between assembly graph edges with weights (if `--debug` was used)
- `<output_dir>/graph_links.tsv` list of graph links between assembly graph edges with weights (if `--debug` was used)

<a name="yaml"></a>
**_Specifying input data with YAML data set file_**

BinSPreader currently supports multiple paired-end or Hi-C libraries described in a YAML file. For example, if you have one paired-end library split into two sets of files

``` bash

    lib_pe1_left_1.fastq
    lib_pe1_right_1.fastq
    lib_pe1_left_2.fastq
    lib_pe1_right_2.fastq
```

and one Hi-C library

``` bash

    lib_hic1_left.fastq
    lib_hic1_right.fastq
```

YAML file should look like this:

``` bash

    [
      {
        orientation: "fr",
        type: "paired-end",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_pe1_right_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_pe1_right_2.fastq" 
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_pe1_left_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_pe1_left_2.fastq" 
        ]
      },
      {
        orientation: "fr",
        type: "paired-end",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_hic1_right.fastq" 
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_hic1_left.fastq"
        ]
      }
    ]
```
