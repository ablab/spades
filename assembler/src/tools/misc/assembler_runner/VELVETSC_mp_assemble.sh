#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


OUTPUT_DIR=RUNNER_PARAM_OUTPUT_DIR

if [[ -z 'RUNNER_PARAM_INTERLACED' ]]; then
  RUNNER_PARAM_VELVET_BASE/shuffleSequences_fastq.pl RUNNER_PARAM_LEFT RUNNER_PARAM_RIGHT $OUTPUT_DIR/shuffled.fastq
  SHUFFLED=$OUTPUT_DIR/shuffled.fastq
else
  SHUFFLED=RUNNER_PARAM_INTERLACED
fi
if [[ -z 'RUNNER_PARAM_MP_INTERLACED' ]]; then
  RUNNER_PARAM_VELVET_BASE/shuffleSequences_fastq.pl RUNNER_PARAM_MP_LEFT RUNNER_PARAM_MP_RIGHT $OUTPUT_DIR/shuffled_mp.fastq
  SHUFFLED_MP=$OUTPUT_DIR/shuffled_mp.fastq
else
  SHUFFLED_MP=RUNNER_PARAM_MP_INTERLACED
fi
RUNNER_PARAM_VELVET_BASE/velveth $OUTPUT_DIR RUNNER_PARAM_K -shortPaired -fastq $SHUFFLED -shortPaired2 -fastq $SHUFFLED_MP
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -ins_length2 RUNNER_PARAM_MP_INSERT_SIZE -ins_length2_sd RUNNER_PARAM_MP_DEVIATION -scaffolding yes 
mv $OUTPUT_DIR/contigs.fa $OUTPUT_DIR/scaffolds.fa
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -ins_length2 RUNNER_PARAM_MP_INSERT_SIZE -ins_length2_sd RUNNER_PARAM_MP_DEVIATION -scaffolding no
