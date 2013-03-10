#!/usr/bin/env bash

REFERENCE_NAME="U00096.2"
REFERENCE_FASTA="U00096.2.fas"
REGION="1-50000"
REFERENCE_SLICE="1-51000"
REFERENCE_REGION_FASTA="ref50k.fas"

## Preparation: #######################################################
sambamba slice ../ecoli_full.bam $REFERENCE_NAME:$REGION -o raw.co.bam
sambamba sort -n --tmpdir=. raw.co.bam -o raw.bam && rm raw.co.bam
bam2fastq raw.bam && mv raw.bam.fastq test.fastq
tmap index -f $REFERENCE_FASTA
fasta_slice $REFERENCE_FASTA $REFERENCE_SLICE > $REFERENCE_REGION_FASTA
tmap index -f $REFERENCE_REGION_FASTA
 
#######################################################################
## Evaluation:

# Run read correction, map corrected reads to the reference, and sorted resulted BAM file by name
hammer-it 2>clusters.txt
tmap mapall -f $REFERENCE_FASTA -r test.fasta -n 16 -o 1 -s corrected.bam -v stage1 map1 map2 map3 map4
sambamba sort -n --tmpdir=. corrected.bam && mv corrected.sorted.bam corrected.bam

# Evaluate clustering
clusters.rb > kmers.fasta
tmap mapall -f $REFERENCE_REGION_FASTA -r kmers.fasta -n 16 -o 1 -s kmers.bam -v stage1 map1 map2 map3 map4
sambamba sort -n --tmpdir=. kmers.bam && mv kmers.sorted.bam kmers.bam
cluster_summaries kmers.bam > cluster_summaries.txt

# Evaluate error correction
comparator raw.bam corrected.bam
sambamba sort --tmpdir=. ruined.bam && mv ruined.sorted.bam ruined.bam
sambamba sort --tmpdir=. damaged.bam && mv damaged.sorted.bam damaged.bam
