#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


K=RUNNER_PARAM_K
CPU=RUNNER_PARAM_THREADS
MIN_OVERLAP=45
ASSEMBLE_OVERLAP=45
MIN_PAIRS=5
RUNNER_PARAM_SGA_BASE/sga preprocess --pe-mode 1 -o reads.pp.fastq frag1 frag2
RUNNER_PARAM_SGA_BASE/sga index --algorithm=ropebwt -t $CPU reads.pp.fastq
RUNNER_PARAM_SGA_BASE/sga correct -k $K -t $CPU -o reads.ec.fastq reads.pp.fastq
RUNNER_PARAM_SGA_BASE/sga index --algorithm=ropebwt -t $CPU reads.ec.fastq
RUNNER_PARAM_SGA_BASE/sga filter -t $CPU reads.ec.fastq
RUNNER_PARAM_SGA_BASE/sga overlap -m $MIN_OVERLAP -t $CPU reads.ec.filter.pass.fa
RUNNER_PARAM_SGA_BASE/sga assemble -o primary reads.ec.filter.pass.asqg.gz
ln -s primary-contigs.fa ctg.fasta
bwa index ctg.fasta
bwa aln -t $CPU ctg.fasta frag1 > frag1.sai
bwa aln -t $CPU ctg.fasta frag2 > frag2.sai
bwa sampe ctg.fasta frag1.sai frag2.sai frag1 frag2 > frag.sam
samtools view -Sb frag.sam > libPE.bam
RUNNER_PARAM_SGA_BASE/sga-bam2de.pl -n $MIN_PAIRS --prefix libPE libPE.bam
RUNNER_PARAM_SGA_BASE/sga-astat.py libPE.bam > libPE.astat
RUNNER_PARAM_SGA_BASE/sga scaffold -m 200 -a libPE.astat -o scf --pe libPE.de ctg.fasta
RUNNER_PARAM_SGA_BASE/sga scaffold2fasta -a primary-graph.asqg.gz -o scf.fasta scf
