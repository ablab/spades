#!/usr/bin/python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing the correct choosing of kmers in SPAdes
# provide a path to .yaml file with test description

import teamcity_workflow


def check_kmer_set(dataset_info, test, output_dir):
    return 0


teamcity_workflow.main(check_test=check_kmer_set)