# SPAligner

Tool for fast and accurate alignment of nucleotide sequences (s.a. long reads, coding sequences, etc.) to assembly graphs. 

## Running SPAligner

    spaligner spaligner_config.yaml \    # config file 
              -d pacbio \                # data type: pacbio, nanopore
              -g assembly_graph.gfa \    # gfa-file with assembly graph 
              -k 77 \                    # graph K-mer size
              -s pacbio_reads.fastq.gz \ # sequences to align in fasta/fastq formats
              -t 8                       # number of threads, 8 by default

By default, spaligner_config.yaml will be installed into /usr/share/spaligner/ or can be found in assembler/projects/spaligner/.

Alignments will be saved to spaligner_result/alignment.tsv by default.


## Compilation

    git clone https://github.com/ablab/spades.git
    cd spades/assembler/
    mkdir build && cd build && cmake ../src
    make spaligner

Now to run SPAligner move to folder `assembler/` and execute

    build/bin/spaligner

## Output

SPAligner can represent the results in three formats: *.tsv (default), *.fasta and [*.gpa](https://github.com/ocxtal/gpa "GPA-format spec").

    spaligner_result/alignment.tsv              tab-separated file with alignments information, each line represents an alignment of a single sequence
    spaligner_result/alignment.fasta            each record represents alignment of a sequence onto assembly graph
    spaligner_result/alignment.gpa              alignment stored in gpa-format


## Results interpretation

SPAligner can represent the results in three formats: *.tsv (default), *.fasta and [*.gpa](https://github.com/ocxtal/gpa "GPA-format spec").
Name of each record in fasta files shows information about alignment position on graph and on sequence. 

**Example 1**

```
>name|Edges=1-,2+,5+|start_g=3283|end_g=35|start_s=0|end_s=291
ATGAAAATCACTCCTGAACAGGCTCGTGAGGCTCTGGATGCCTGGATATGTCGACCAGGAATGACACAGGAGCAGGCGACGATATTAATCACTGAAGCATTCTGGGCTTTGAAAGAGCGCCCGAACATCGATGTTCAGCGTGTCACATATGAAGGTGGCGCGATTGATCAGCGAGCGCTTGGCGTTAATCGAGTGAAGATATTTGAACGCTGGAAGGCTATCGACACCAGGGATAAGC
GTGAAAAGTTCACGGCGCTAGTGCCTGCAATTATGGAGGCTACCACTGGATGA
```

name — sequence name<br/>
1-,2+,5+ — alignment path<br/>
3283 — start position of alignment on the first edge of the path (here on conjugate edge to edge with id=1)<br/>
35 —  end position of alignment on the last edge of the path (here on edge with id=5)<br/>
0 — start position of alignment on sequence<br/>
291 — end position of  alignment on sequence<br/>


Each line in tsv-file represents alignments of a single read.

**Example 2**

```
name    0      2491  536  1142  2491	44+,24+,22+,1+,38-	909,4,115,1,1142	AAACTTTTATTGTGCATACGGCGATTAAGACGGGAAAAGTCGGTGAT...
```

name — sequence name<br/>
0 — start position of alignment on sequence<br/>
2491 — end position of  alignment on sequence<br/>
536 — start position of alignment on the first edge of the Path (here on edge with id=44)<br/>
1142 —  end position of alignment on the last edge of the Path (here on conjugate edge to edge with id=38)<br/>
2491 — sequence length<br/>
44+,24+,22+,1+,38- — Path of the alignment <br/>
909,4,115,1,1142 — lengths of the alignment on each edge of the Path respectively (44+,24+,22+,1+,38-) <br/>
AGGTTGTTTTTTGTTTCTTCCGC... — sequence of alignment Path <br/>



**Example 3**<br/>

Sometimes sequence alignment on the graph can be represented as several non-overlapping subpaths (if there is no alignment with appropriate score between two consecutive bwa anchors). <br/>
So, there can be several unconnected alignments of sequence onto assembly graph and several start positions, end positions, paths etc.:

```
name     4,10  7,19       2,7   5,6       19    123+;288-,128+       3;3,6 GAT;TTATCCGGG
```

The sequence *name* has two alignments on the graph:

1. The first alignment starts on sequence on position 4 and ends on position 7.
Corresponding path consists of a single edge 123+ (i.e. 123+) with start on position 2 and end on position 5.
Path sequence: GAT.
2. While the second alignment covers the end of the sequence, starting on sequence on position 10 and ending on position 19. 
Corresponding path consists of two edges 288- and 128+ (i.e. 288-,128+) with start on position 7 (on 288-) and end on position 6 (on 128+).
Path sequence: TTATCCGGG.

If a sequence was not fully aligned, SPAligner tries to prolong the longest alignment subpath in order to reconstruct a full alignment path. In **Example 2** SPAligner was not able to prolong any of two given alignments.


## Algorithm overview

Overview of the alignment of the nucleotide query sequence *S* (orange bar) to assembly graph *G*. Assembly graph edges are considered directed left-to-right (explicit edge orientation was omitted to improve the clarity).

![pipeline](pipeline.jpg)

1. **Anchor search.** Anchors (regions of high similarity) between the query and the edge labels are identified with [BWA-MEM](http://bio-bwa.sourceforge.net/). 
2. **Anchor filtering.** Anchors shorter than *K*, assembly graph *K*-mer size,(anchors 2, 6, 11), anchors “in the middle” of long edge (anchor 7) or ambiguous anchors (anchor 10 mostly covered by anchor 9, both anchors 4 and 5) are discarded.
3. **Anchor chaining.** Heaviest chain of compatible anchors (chain 3->8->9) is determined.
4. **Reconstruction of filling paths.** Paths for fragments of the query between the consecutive chain anchors (as well as left- and right-most fragments) are reconstructed. The procedure is performed using fast library for sequence alignment [Edlib](https://github.com/Martinsos/edlib).

## Parameters tuning

Full list of parameters can be found in spaligner_config.yaml. 

### General alignment parameters

* `run_dijkstra: true` Run Dijkstra algorithm to find alignment between anchors, if `run_dijkstra: false`, SPAligner will check limited number of paths and return the best one.
* `restore_ends: true` Restore alignment path before leftmost anchor and after rightmost anchor.


* `internal_length_cutoff: 200` Anchors with length < `internal_length_cutoff` will be filtered out.
* `path_limit_stretching: 1.3` Pair of anchors is considered to be compatible if (the minimal distance between them in graph) 
                                < `path_limit_stretching` * (the distance between their positions on sequence).
* `path_limit_pressing: 0.7` Pair of anchors is considered to be compatible if (the minimal distance between them in graph) 
                                > `path_limit_pressing` * (the distance between their positions on sequence).
* `max_path_in_chaining: 15000` Limit on number of paths to consider between two anchors on anchors chaining step.
* `max_vertex_in_chaining: 5000` Limit on number of vertices to consider between two anchors on anchors chaining step.

### Internal run parameters

* `queue_limit: 1000000` Limit on queue length. 
* `iteration_limit: 1000000` Limit on total number of queue extraction. 
* `updates_limit: 1000000` Limit on number of updates of shortest distance for all states.
* `find_shortest_path: true` If `find_shortest_path: false` Dijkstra algorithm will stop when it reaches finish state without searching for shortest path.
* `restore_mapping: false` If `restore_mapping: true` Dijkstra algorithm will return full alignment information (used by for developers).
* `penalty_ratio: 0.1` Algorithm never considers states representing an alignment of a prefix S[0:i] with score more than min_score(i) + i*penalty_ratio (for nucleotide sequence alignment only).
* `max_ed_proportion: 3` Maximal edit distance is bounded by a fraction of the query sequence length |S|/max_ed_proportion. Increase of max_ed_proportion leads to shorter alignments but with higer identity.
* `ed_lower_bound: 500` Minimal penalty score of alignment.
* `ed_upper_bound: 2000` Maximal penalty score of alignment.
* `max_gs_states: 120000000` If number of queue states exceeds `max_gs_limit` then shortest path search is not performed (for nucleotide sequence alignment only).
* `max_restorable_length: 5000` If distance between two anchors or between leftmost/rightmost anchor and start/end exceeds `max_restorable_length` then shortest path search is not performed.

Increase of `max_gs_states`, `max_restorable_length`, `queue_limit`, `iteration_limit` or `updates_limit` may lead to longer alignments with the same identity level, but slows down the process and can use much more memory. Please change them if you 100% confident in what you are doing.

Turning off restore_ends or run_dijkstra in nucleotide sequence alignment mode leads to shorter alignments, but considerable speed-up.


## Contacts

For any questions or suggestions please do not hesitate to contact Tatiana Dvorkina <tedvorkina@gmail.com>.
