# SPAligner

Tool for fast and accurate alignment of nucleotide sequences (s.a. long reads, coding sequences, etc.) to assembly graphs. 


## Compilation
---
    git clone https://github.com/ablab/spades.git
    cd algorithmic-biology/assembler/
    ./prepare_cfg
    make -C build/release/projects/spaligner/ -j4 longreads_aligner


## Running SPAligner
---

To align PacBio reads realpb.fasta (accepts fasta/fastq files) to *E.coli* assembly graph built for K=77 (dataset is available [here](https://figshare.com/s/004baf22fc1bfd758f5b "Figshare DB")):
``` 
./build/release/bin/longreads_aligner ./src/projects/spaligner/config.yaml -K 77 -d pacbio -g ecoli.gfa -s realpb.fasta -o test_ecoli > test_ecoli.log
```
Alignments will be saved to test_ecoli.tsv. 

Run *help* to see full list of options:
```
./build/release/bin/longreads_aligner -h
```


## Results interpretation
---
SPAligner can represent the results in three formats: *.tsv (default), *.fasta and [*.gpa](https://github.com/ocxtal/gpa "GPA-format spec").
Each line in tsv-file represents alignments of a single read.

**Examples:**

* Example1

```
name    0      24911  53661  11429  24911	444-,30+,494-,264-,342-,264-,342-,264-,342-,264-,361-,224+,225+,1+,386-	9096,23,157,41,140,41,140,41,140,41,1573,4,1156,1,11429	AAACTTTTATTGTGCATACGGCGATTAAGACGGGAAAAGTCGGTGAT...
```

Description:
name -- read name
0 -- start position of alignment on read
24911 -- end position of alignment on read
53661 -- start position of alignment on the first edge of the Path (here on conjugate edge to edge with id=444)
11429 --  end position of alignment on the last edge of the Path (here on conjugate edge to edge with id=386)
24911 -- read length
444-,30+,494-,264-,342-,264-,342-,264-,342-,264-,361-,224+,225+,1+,386- -- Path of the alignment
9096,23,157,41,140,41,140,41,140,41,1573,4,1156,1,11429 -- lengths of the alignment on each edge of the Path respectively (444-,30+,494-,264-,342-,...)
AGGTTGTTTTTTGTTTCTTCCGC... -- sequence of alignment Path


* Example2
Sometimes read alignment on the graph can be represented as several non-overlapping subpaths (if there is no alignment with appropriate score between two consecutive bwa hits).
So, there can be several unconnected alignments of read onto assembly graph and several start positions, end positions, paths etc.:

```
name2     4,10  7,19       2,7   5,6       19    123+;288-,128+       3;3,6 GAT;TTATCCGGG
```

The read *name2* has two alignments on the graph:

1. The first alignment starts on read on position 4 and ends on position 7.
Corresponding path consists of a single edge 123+ (i.e. 123+) with start on position 2 and end on position 5.
Path sequence: GAT
2. While the second alignment covers the end of the read, starting on read on position 10 and ending on position 19. 
Corresponding path consists of two edges 288- and 128+ (i.e. 288-,128+) with start on position 7 (on 288-) and end on position 6 (on 128+).
Path sequence: TTATCCGGG


If read was not fully aligned, SPAligner tries to prolong the longest alignment subpath in order to reconstruct full alignment path. In *Example2* SPAligner was not ably to prolong any of two given alignments.



## Future plans 
---
1. Add amino acid sequence support.
2. Alignment speed-up.
