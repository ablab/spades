#!/bin/sh

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


RUNNER_PARAM_CABOG_BASE/fastqToCA -insertsize RUNNER_PARAM_INSERT_SIZE RUNNER_PARAM_DEVIATION -type 'RUNNER_PARAM_PHRED_TYPE' -libraryname reads -mates RUNNER_PARAM_LEFT,RUNNER_PARAM_RIGHT >reads.frg 2>>RUNNER_PARAM_STDERR_LOG
RUNNER_PARAM_CABOG_BASE/runCA -d . -p asm -s RUNNER_PARAM_CABOG_CONFIG reads.frg >> RUNNER_PARAM_STDOUT_LOG 2>> RUNNER_PARAM_STDERR_LOG

 
