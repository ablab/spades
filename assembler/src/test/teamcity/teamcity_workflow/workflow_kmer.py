#!/usr/bin/python

############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing the correct choosing of kmers in SPAdes
# provide a path to .yaml file with test description

import workflow_base


def check_kmer_set(dataset_info, test, output_dir):
    return 0


workflow_base.main(check_test=check_kmer_set)