#!/usr/bin/python

############################################################################
# Copyright (c) 2020 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing the correct choosing of kmers in SPAdes
# provide a path to .yaml file with test description

import workflow_base
from workflow_base import log
import os
import sys
from site import addsitedir

spades_home = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))), "..", "..", "..", "..")
ext_python_modules_home = os.path.join(spades_home, "ext", "src", "python_libs")
addsitedir(ext_python_modules_home)
if sys.version.startswith("2."):
    import pyyaml2 as pyyaml
elif sys.version.startswith("3."):
    import pyyaml3 as pyyaml


def get_kmer_list(path):
    run_spades_yaml = os.path.join(path, "out", "run_spades.yaml")
    stages = pyyaml.load(open(run_spades_yaml))
    kmers = []
    for stage in stages:
        stage_name = stage['STAGE']
        if stage_name[0] == 'K' and stage_name[1:].isdigit():
            kmers.append(stage_name)
    return kmers


def check_kmer_set(dataset_info, test, output_dir):
    if 'etalon_saves' in dataset_info:
        log.log("Checking the correctness of K")
        etalon_folder = dataset_info["etalon_saves"]
        if ("name" in test):
            etalon_folder += test["name"]

        out_klist = get_kmer_list(output_dir)
        log.log("Output kmers: " + str(out_klist))

        etalon_klist = get_kmer_list(os.path.join(etalon_folder))
        log.log("Etalon kmers: " + str(etalon_klist))

        if out_klist != etalon_klist:
            log.err("Kmers in output(" + str(out_klist) + ") != kmers in etalon(" + str(etalon_klist) + ")")
            return 12
        return 0
    else:
        log.err("Etalon folder wasn't set in test config!")
        return 12

workflow_base.main(check_test=check_kmer_set)