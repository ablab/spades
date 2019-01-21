# SPAligner

Tool for fast and accurate alignment of nucleotide sequences (s.a. long reads, coding sequences, etc.) to assembly graphs. 


## Compilation
---
    git clone https://github.com/ablab/algorithmic-biology.git
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
SPAligner can represent the results in two formats: *.tsv (default) and [*.gpa](https://github.com/ocxtal/gpa "GPA-format spec").
Each line in tsv-file represents alignments of a single read.

**Examples:**

* Example1

```
name    0,      24911,  53661,  11429,  24911   444-,30+,494-,264-,342-,264-,342-,264-,342-,264-,361-,224+,225+,1+,386-,;       9096,23,157,41,140,41,140,41,140,41,1573,4,1156,1,11429,;       444- (53705,62757) [47,9488], 361- (0,1573) [10253,11907], 225+ (0,1156) [11912,13106], 386- (0,11330) [13106,24812], ; AAACTTTTATTGTGCATACGGCGATTAAGACGGGAAAAGTCGGTGAT...
```

Description:
name -- Read name
0, -- start position of alignment on Read
24911, -- end position of alignment on Read
53661, -- start position of alignment on the first edge of the Path (here on conjugate edge to edge with id=444)
11429, --  end position of alignment on the last edge of the Path (here on conjugate edge to edge with id=386)
24911 -- Read length
444-,30+,494-,264-,342-,264-,342-,264-,342-,264-,361-,224+,225+,1+,386-,; -- Path of the alignment
9096,23,157,41,140,41,140,41,140,41,1573,4,1156,1,11429,; -- lengths of the alignment on each edge of the Path respectively (444-,30+,494-,264-,342-,...)
444- (53705,62757) [47,9488], 361- (0,1573) [10253,11907], 225+ (0,1156) [11912,13106], 386- (0,11330) [13106,24812], ; -- BWA-MEM hits on edges
AGGTTGTTTTTTGTTTCTTCCGC... -- sequence of alignment Path



* Example2

As there can be several alignments of Read onto assembly graph, there can be several start positions, end positions, paths etc.:

```
name2     4,0,  17,19,       2,44,   15,10,       19    123+,;288-,128+,;       13,;10,10,; 123+ (2,15) [4,17], ;128+ (3,10) [12,19], ;  
GATACGGTTATCC;GGGCGATACGGTTATCCGGG
```

The read *name2* has two alignments on the graph:

1. The first alignment starts on read on position 4 and ends on position 17. 
Corresponding path consists of a single edge 123+ (i.e. 123+,;) with start on position 2 and end on position 15.
Path sequence: GATACGGTTATCC
2. While the second alignment is a "full" alignment: from position 0 to position 19 on the read.
Corresponding path consists of two edges 288- and 128+ (i.e. 288-,128+;) with start on position 44 and end on position 10.
Path sequence: GGGCGATACGGTTATCCGGG


## Supported assembly graph types
---
SPAligner accepts graphs in *.gfa format, interpreting `S` segments as edges of the condensed de Bruijn graph (and `L` links -- as vertices).


## Future plans 
---
1. Add amino acid sequence support.
2. Alignment speed-up.
