# SPAdes output

SPAdes stores all output files in `<output_dir> `, which is set by the user.

-   `<output_dir>/corrected/` directory contains reads corrected by BayesHammer in `*.fastq.gz` files; if compression is disabled, reads are stored in uncompressed `*.fastq` files
-   `<output_dir>/scaffolds.fasta` contains resulting scaffolds (recommended for use as resulting sequences)
-   `<output_dir>/contigs.fasta` contains resulting contigs
-   `<output_dir>/assembly_graph_with_scaffolds.gfa` contains SPAdes assembly graph and scaffolds paths in [GFA 1.2 format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md)
-   `<output_dir>/assembly_graph.fastg` contains SPAdes assembly graph in [FASTG format](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf)
-   `<output_dir>/contigs.paths` contains paths in the assembly graph corresponding to contigs.fasta (see details below)
-   `<output_dir>/scaffolds.paths` contains paths in the assembly graph corresponding to scaffolds.fasta (see details below)

## Contigs and scaffolds format

Contigs/scaffolds names in SPAdes output FASTA files have the following format:
`>NODE_3_length_237403_cov_243.207`
Here `3` is the number of the contig/scaffold, `237403` is the sequence length in nucleotides and `243.207` is the k-mer coverage for the last (largest) k value used. Note that the k-mer coverage is always lower than the read (per-base) coverage.

In general, SPAdes uses two techniques for joining contigs into scaffolds. First one relies on read pairs and tries to estimate the size of the gap separating contigs. The second one relies on the assembly graph: e.g. if two contigs are separated by a complex tandem repeat that cannot be resolved exactly, contigs are joined into a scaffold with a fixed gap size of 100 bp. Contigs produced by SPAdes do not contain N symbols.

## Assembly graph formats
SPAdes produces assembly graph in GFA 1.2 and legacy FASTG formats. 
To view GFA and FASTG files we recommend to use [Bandage-NG visualization tool](https://github.com/asl/BandageNG).

### GFA
SPAdes encodes contigs before repeat resolution as [segments](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#s-segment-line) in the GFA assembly graph. Their overlaps and corresponding connections are represented as [links](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#l-link-line) in the graph. Each segment is complemented with two tags: `DP`, which represents the average k-mer coverage depths, and `KC`, which encodes essentially the same value but in the raw number of k-mers from which this segment was assembled. Both coverages are calculated using the last used k-mer length.

In addition to segments and links, SPAdes encodes scaffold paths through the assembly graph using GFA [path records](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#p-path-line). The path name is same as a scaffold name. Also, if the scaffold is circular, the `TP:Z:circular` tag is added.

Note that scaffolds generally represent gapped paths through the assembly graph. To be able to represent this, a GFA v1.2 feature called [jump links](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#extension-to-use-jump-connections-since-v12) is employed.

If the consumer of an assembly graph is not ready to cope with the GFA v1.2 format, it is possible to switch to legacy GFA 1.1 format used in previous SPAdes versions using the `--gfa11` option. In GFA v1.1 / v1.0 mode it is not possible to represent scaffolds paths fully as they must be contiguous. To represent scaffolds gaps in GFA v1.0-style assembly graph, the scaffold is broken down into contigs and each contiguous segment is represented as a separate GFA path. The name of the GFA path includes the original scaffold name as well as a segment number.

To further showcase the differences between GFA v1.2 and v1.1 / v1.0-style GFA paths, here is the snippet from the real output graphs:

``` plain
   # GFA v1.2
   J    16337   +    13759    -   *   SC:i:1   
   P    NODE_16_length_112493_cov_221.698741      16337+;13759-,4559-,22173-    *

   # GFA v1.1 / 1.0
   P    NODE_16_length_112493_cov_221.698741_1    16337+    *
   P    NODE_16_length_112493_cov_221.698741_2    13759-,4559-,22173-   *
```

Note that the `_1` name suffix is always present in GFA v1.1 / 1.0 version even for the scaffolds that do not have any gaps inside.

### FASTG
The sequences stored in `assembly_graph.fastg` correspond to contigs before repeat resolution (edges of the assembly graph). Paths corresponding to contigs after repeat resolution (scaffolding) are stored in `contigs.paths` (`scaffolds.paths`) in the format accepted by Bandage (see [Bandage wiki](https://github.com/rrwick/Bandage/wiki/Graph-paths) for details). The example is given below.

Let the contig with the name `NODE_5_length_100000_cov_215.651` consist of the following edges of the assembly graph:

``` plain
    >EDGE_2_length_33280_cov_199.702
    >EDGE_5_length_84_cov_321.414"
    >EDGE_3_length_111_cov_175.304
    >EDGE_5_length_84_cov_321.414"
    >EDGE_4_length_66661_cov_223.548
```

Then, `contigs.paths` will contain the following record:

``` plain
    NODE_5_length_100000_cov_215.651
    2+,5-,3+,5-,4+
```


Since the current version of Bandage does not accept paths with gaps, paths corresponding contigs/scaffolds jumping over a gap in the assembly graph are split by semicolon at the gap positions. For example, the following record

``` plain
    NODE_3_length_237403_cov_243.207
    21-,17-,15+,17-,16+;
    31+,23-,22+,23-,4-
```

states that `NODE_3_length_237403_cov_243.207` corresponds to the path with 10 edges, but jumps over a gap between edges `EDGE_16_length_21503_cov_482.709` and `EDGE_31_length_140767_cov_220.239`.

## Complete list of output files

The full list of `<output_dir>` content is presented below:

- scaffolds.fasta - resulting scaffolds (recommended for use as resulting sequences)
- contigs.fasta - resulting contigs
- assembly_graph.fastg - assembly graph
- contigs.paths - contigs paths in the assembly graph
- scaffolds.paths - scaffolds paths in the assembly graph
- before_rr.fasta - contigs before repeat resolution

- corrected/ - files from read error correction
    - configs/ - configuration files for read error correction
    - corrected.yaml - internal configuration file
    - Output files with corrected reads

- params.txt - information about SPAdes parameters in this run
- spades.log - SPAdes log
- dataset.info - internal configuration file
- input_dataset.yaml - internal YAML data set file
- K<##>/ - directory containing intermediate files from the run with K=<##>. These files should not be used as assembly results; use resulting contigs/scaffolds in files mentioned above.


SPAdes will overwrite these files and directories if they exist in the specified `<output_dir>`.

## plasmidSPAdes output

plasmidSPAdes and metaplasmidSPAdes output only DNA sequences from putative plasmids. Output file names and formats remain the same as in SPAdes (see [previous](output.md#spades-output) section), with the following differences.

For all plasmidSPAdes' contig names in `contigs.fasta`, `scaffolds.fasta` and `assembly_graph.fastg` we append suffix `_component_X`, where `X` is the id of the putative plasmid, which the contig belongs to. Note that plasmidSPAdes may not be able to separate similar plasmids and thus their contigs may appear with the same id.


## metaplasmidSPAdes and metaviralSPAdes output
The repeat resolution and extrachromosomal element detection in metaplasmidSPAdes/metaviralSPAdes is run independently for different coverage cutoffs values (see [paper](https://genome.cshlp.org/content/29/6/961.short) for details). In order to distinguish contigs with putative plasmids detected at different cutoff levels we extend the contig name in FASTA file with cutoff value used for this particular contig (in format `_cutoff_N`). This is why, in the contrast to regular SPAdes pipeline, there might be a contig with `NODE_1_` prefix for each cutoff with potential plasmids detected. In the following example, there were detected two potential viruses using cutoff 0, one virus was detected with cutoff 5 and one with cutoff 10. We also add a suffix that shows the structure of the suspective extrachromosomal element.

In the metaplasmid mode SPAdes outputs only circular putative plasmids.
In the metaviral mode SPAdes also outputs linear putative viruses and linear viruses with simple repeats ('9'-shaped components in the assembly graph) sequences.

``` 
>NODE_1_length_40003_cov_13.48_cutoff_0_type_circular
>NODE_2_length_30000_cov_4.20_cutoff_0_type_linear
>NODE_1_length_20000_cov_20.42_cutoff_5_type_circular
>NODE_1_length_10000_cov_198.4_cutoff_10_type_linearrepeat
```

## biosyntheticSPAdes output

biosyntheticSPAdes outputs four files of interest:
- scaffolds.fasta – contains DNA sequences from putative biosynthetic gene clusters (BGC). Since each sample may contain multiple BGCs and biosyntheticSPAdes can output several putative DNA sequences for each cluster, for each contig name we append suffix `_cluster_X_candidate_Y`, where X is the id of the BGC and Y is the id of the candidate from the BGC.
- raw_scaffolds.fasta - SPAdes scaffolds generated without domain-graph related algorithms. Very close to the regular scaffolds.fasta file.
- hmm_statistics.txt - contains statistics about BGC composition in the sample. First, it outputs the number of domain hits in the sample. Then, for each BGC candidate we output domain order with positions on the corresponding DNA sequence from scaffolds.fasta.
- domain_graph.dot - contains domain graph structure that can be used to assess complexity of the sample and structure of BGCs. For more information about domain graph construction, please refer to the [biosyntheticSPAdes](https://genome.cshlp.org/content/early/2019/06/03/gr.243477.118?top=1) paper.


## rnaSPades output

See [rnaSPAdes section](rna.md#rnaspades-output).


## Genome assembly evaluation

[QUAST](https://quast.sourceforge.net/) may be used to generate summary statistics (N50, maximum contig length, GC %, \# genes found in a reference list or with built-in gene finding tools, etc.) for a single assembly. It may also be used to compare statistics for multiple assemblies of the same data set (e.g., SPAdes run with different parameters, or several different assemblers).

