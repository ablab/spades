# BinSPreader - binning refining using assembly graphs

## MANUAL

### Description

BinSPreader is a novel tool that attempts to refine metagenome-assembled genomes (MAGs) obtained from existing tools. BinSPreader exploits the assembly graph topology and other connectivity information, such as paired-end and Hi-C reads, to refine the existing binning, correct binning errors, propagate binning from longer contigs to shorter contigs and infer contigs belonging to multiple bins.

### Installation 

``` 
cd spades/assembler/
mkdir build && cd build && cmake ../src
make bin-refine
```
Now to run BinSPreader move to folder `assembler/` and execute 

`build/bin/hicspades-binner`

### Input

The tool has two mandatory options: 
- Assembly graph file in [GFA 1.0 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), with scaffolds included as path lines. Alternatively, scaffolds can be provided separately using `--path` option. 
-Initial 

Synopsis: `hicspades-binner <graph (in GFA)> <dataset description (in YAML)> <output directory> [OPTION...]`

The options are:

`-t, --threads <int> `
    # of threads to use

`-e, --enzymes <string> `
    Comma-separated string of restriction enzyme recognition sites

`--tmp-dir <dir name> `
    scratch directory to use

`--min-ctg-len <int> `
    Minimum contig length for binning

`--path-links-thr <int> `
    Minimum total number of links between contigs

`--edge-links-thr <int>` 
    Minimum number of links between long edges

`-h, --help `
    print help message

<a name="yaml"></a>
**_Specifying input data with YAML data set file_**

hicSPAdes-binner currently supports a single Hi-C library described in a YAML file. For example, if your Hi-C library is split into two pairs of files

``` bash

    lib_hic_left_1.fastq
    lib_hic_right_1.fastq
    lib_hic_left_2.fastq
    lib_hic_right_2.fastq
```

YAML file should look like this:

``` bash

    [
      {
        orientation: "fr",
        type: "hic",
        right reads: [
          "/FULL_PATH_TO_DATASET/lib_hic_right_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_hic_right_2.fastq" 
        ],
        left reads: [
          "/FULL_PATH_TO_DATASET/lib_hic_left_1.fastq",
          "/FULL_PATH_TO_DATASET/lib_hic_left_2.fastq" 
        ]
      }
    ]
```

### Output

hicSPAdes-binner stores all output files in `<output_dir> `, which is set by the user.

-   `<output_dir>/clustering.mcl` contains resulting scaffold clustering in MCL format
-   `<output_dir>/clustering.tsv` contains resulting scaffold clustering in TSV format
-   `<output_dir>/basic_stats.tsv` contains various per-cluster statistics
-   `<output_dir>/contact_map.tsv` contains hicSPAdes scores between input scaffolds, as well as other scaffold statistics
