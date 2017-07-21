#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import logging
import os
import shutil
import sys
from os.path import abspath, dirname, realpath, join, isfile
from site import addsitedir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result_scaffolds_filename",
                    help="path to file with result scaffolds",
                    action="store")
    parser.add_argument("--assembled_scaffolds_filename",
                    help="path to file with assembled scaffolds",
                    action="store")
    parser.add_argument("--bin_home",
                    help="path to bin home",
                    action="store")
    parser.add_argument("--ext_python_modules_home",
                    help="path to ext python modules home",
                    action="store")
    parser.add_argument("--output_dir",
                    help="path to output dir",
                    action="store")
    parser.add_argument("--truseq_long_reads_file_base",
                    help="path to file with truseq long reads",
                    action="store")
    parser.add_argument("--dataset_yaml_file",
                    help="path to yaml file with dataset",
                    action="store")
    parser.add_argument("--threads",
                    type=int,
                    help="number of threads",
                    action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    # create logger
    log = logging.getLogger("Postprocessing")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    addsitedir(args.ext_python_modules_home)
    # save dataset from yaml
    if sys.version.startswith("2."):
        import pyyaml2 as pyyaml
    elif sys.version.startswith("3."):
        import pyyaml3 as pyyaml

    dataset_data = pyyaml.load(open(args.dataset_yaml_file))

    # init python_lib_folder
    python_modules_home = abspath(dirname(realpath(__file__)))
    source_dirs = ["..", "../truspades", "../common", "../executors"]
    for dir_name in source_dirs:
        sys.path.append(join(python_modules_home, dir_name))

    # import alignment and molecule_postprocassing
    import alignment
    import moleculo_postprocessing
    #  run command
    if isfile(args.result_scaffolds_filename):
        shutil.move(args.result_scaffolds_filename, args.assembled_scaffolds_filename)
    alignment_bin = os.path.join(args.bin_home, "spades-bwa")
    alignment_dir = os.path.join(args.output_dir, "alignment")
    sam_files = alignment.align_bwa(alignment_bin, args.assembled_scaffolds_filename,
                                    dataset_data, alignment_dir, log, args.threads)

    moleculo_postprocessing.moleculo_postprocessing(args.assembled_scaffolds_filename,
                                                    args.truseq_long_reads_file_base,
                                                    sam_files, log)


if __name__ == "__main__":
    main()
