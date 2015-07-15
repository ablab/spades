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
RUNNER_PARAM_VELVET_BASE/velveth $OUTPUT_DIR RUNNER_PARAM_K -shortPaired -fastq $SHUFFLED
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -scaffolding yes
mv $OUTPUT_DIR/contigs.fa $OUTPUT_DIR/scaffolds.fa
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -scaffolding no
