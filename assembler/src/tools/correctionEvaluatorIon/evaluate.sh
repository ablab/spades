#!/usr/bin/env bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


REFERENCE_NAME="U00096.2"
REFERENCE_FASTA="U00096.2.fas"
REGION="1-50000"
REFERENCE_SLICE="1-51000"
REFERENCE_REGION_FASTA="ref50k.fas"

## Preparation: #######################################################
sambamba slice ../ecoli_full.bam $REFERENCE_NAME:$REGION -o raw.co.bam
sambamba sort -n -t8 --tmpdir=. raw.co.bam -o raw.bam && rm raw.co.bam
bam2fastq raw.bam && mv raw.bam.fastq test.fastq
tmap index -f $REFERENCE_FASTA
fasta_slice $REFERENCE_FASTA $REFERENCE_SLICE > $REFERENCE_REGION_FASTA
tmap index -f $REFERENCE_REGION_FASTA
 
#######################################################################
## Evaluation:

# Run read correction, map corrected reads to the reference, and sorted resulted BAM file by name
hammer-it hammer-it.cfg
tmap mapall -f $REFERENCE_FASTA -r test.fasta -n 16 -o 1 -s corrected.bam -v stage1 map1 map2 map3 map4
sambamba sort -n -t8 --tmpdir=. corrected.bam && mv corrected.sorted.bam corrected.bam

# Evaluate clustering
tmap mapall -f $REFERENCE_REGION_FASTA -r kmers.fasta -n 16 -o 1 -s kmers.bam -v stage1 map1 map2 map3 map4
sambamba sort -n -t8 --tmpdir=. kmers.bam && mv kmers.sorted.bam kmers.bam
cluster_summaries kmers.bam > cluster_summaries.txt

# Evaluate error correction
comparator raw.bam corrected.bam
sambamba sort -t8 --tmpdir=. ruined.bam && mv ruined.sorted.bam ruined.bam
sambamba sort -t8 --tmpdir=. damaged.bam && mv damaged.sorted.bam damaged.bam
