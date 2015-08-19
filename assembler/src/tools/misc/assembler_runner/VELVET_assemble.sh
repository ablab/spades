#!/bin/bash

OUTPUT_DIR=RUNNER_PARAM_OUTPUT_DIR

RUNNER_PARAM_VELVET_BASE/velveth $OUTPUT_DIR RUNNER_PARAM_K -shortPaired -separate -fastq RUNNER_PARAM_LEFT RUNNER_PARAM_RIGHT
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -scaffolding yes
mv $OUTPUT_DIR/contigs.fa $OUTPUT_DIR/scaffolds.fa
RUNNER_PARAM_VELVET_BASE/velvetg $OUTPUT_DIR -exp_cov auto -cov_cutoff auto -ins_length RUNNER_PARAM_INSERT_SIZE -ins_length_sd RUNNER_PARAM_DEVIATION -scaffolding no
