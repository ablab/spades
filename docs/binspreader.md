# Binning refining using assembly graphs

BinSPreader is a tool that attempts to refine metagenome-assembled genomes
(MAGs) obtained from existing tools. BinSPreader exploits the assembly graph
topology and other connectivity information, such as paired-end and Hi-C reads,
to refine the existing binning, correct binning errors, and propagate binning from
longer contigs to shorter contigs, and infer contigs belonging to multiple MAGs. 
Please refer to the [BinSPreader paper](https://www.sciencedirect.com/science/article/pii/S2589004222010422)
for more details. In addition to increasing the completeness of the bins, refinement 
also enriches bins with contigs containing important conservative genes using the 
short assembly graph edges which are typically underrepresented in state-of-the-art
contig binning methods.

The tool requires initial binning to refine, as well as an assembly graph as a
source of information for refining. Optionally, BinSPreader can be provided with
multiple Hi-C and/or paired-end libraries. The [BinSPreader protocol](https://star-protocols.cell.com/protocols/2802) contains more detailed
instructions on installing and running BinSPreader.

## Compilation

To compile BinSPreader, run

```
./spades_compile -SPADES_ENABLE_PROJECTS=binspreader
```

After the compilation is complete, `binspreader` executable will be located in the `bin` folder.

## Command line options

Required positional arguments: 

- Assembly graph file in [GFA 1.0
  format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md), with
  scaffolds included as path lines. Alternatively, scaffold paths can be
  provided separately using `--path` option in the `.paths` format accepted by
  Bandage (see [Bandage
  wiki](https://github.com/rrwick/Bandage/wiki/Graph-paths) for details).
- Binning output from an existing tool (in `.tsv` format)

### Synopsis
```bash
binspreader <graph (in GFA)> <binning (in .tsv)> <output directory> [OPTION...]
```

### Main options

`--paths`
    provide contigs paths from file separately from GFA

`--dataset` 
    Dataset in [YAML format](running.md#specifying-multiple-libraries-with-yaml-data-set-file) describing Hi-C and paired-end reads

 `-t` 
    Number of threads to use (default: 1/2 of available threads)

 `-m` 
    Allow multiple bin assignment (default: false)
    
 `-Smax|-Smle` 
     Simple maximum or maximum likelihood binning assignment strategy (default: max likelihood)
     
 `-Rcorr|-Rprop` 
     Select propagation or correction mode (default: correction)
     
`--cami` 
    Use CAMI bioboxes binning format
    
`--zero-bin` 
    Emit zero bin for unbinned sequences
    
`--tall-multi` 
    Use tall table for multiple binning result
    
`--bin-dist` 
    Estimate pairwise bin distance (could be slow on large graphs!)
    
`-la` 
    Labels correction regularization parameter for labeled data (default: 0.6)


## Output
BinSPreader stores all output files in the output directory `<output_dir> ` set by the user.

- `<output_dir>/binning.tsv` contains refined binning in `.tsv` format
- `<output_dir>/bin_stats.tsv` contains various per-bin statistics
- `<output_dir>/bin_weights.tsv` contains soft bin weights per contig
- `<output_dir>/edge_weights.tsv` contains soft bin weights per edge

In addition

- `<output_dir>/bin_dist.tsv` contains refined bin distance matrix (if `--bin-dist` was used)
- `<output_dir>/bin_label_1.fastq, <output_dir>/bin_label_2.fastq` read set for bin labeled by `bin_label` (if `--reads` was used)
- `<output_dir>/pe_links.tsv` list of paired-end links between assembly graph edges with weights (if `--debug` was used)
- `<output_dir>/graph_links.tsv` list of graph links between assembly graph edges with weights (if `--debug` was used)


## References

If you are using **BinSPreader** in your research, please cite:

[Tolstoganov et al., 2022](https://www.cell.com/iscience/pdf/S2589-0042(22)01042-2.pdf) and
[Ochkalova et al., 2023](https://www.sciencedirect.com/science/article/pii/S2666166723003842). 