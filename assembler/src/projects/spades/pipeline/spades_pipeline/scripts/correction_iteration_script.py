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
from os.path import abspath, dirname, realpath, join

python_modules_home = abspath(dirname(realpath(__file__)))
sys.path.append(join(python_modules_home, ".."))
import support


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corrected",
                    help="path to file with corrected contigs/scaffolds",
                    action="store")
    parser.add_argument("--assembled",
                    help="path to file with assembled contigs/scaffolds",
                    action="store")
    parser.add_argument("--assembly_type",
                    choices=["contigs", "scaffolds"],
                    help="assembly type: contigs/scaffolds",
                    action="store")
    parser.add_argument("--output_dir",
                    help="path to output dir",
                    action="store")
    parser.add_argument("--bin_home",
                    help="path to bin home",
                    action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    # create logger
    log = logging.getLogger("Mismatch correction " + args.assembly_type)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    # moving assembled contigs (scaffolds) to misc dir
    if os.path.isfile(args.corrected):
        shutil.move(args.corrected, args.assembled)

    # TODO can check only here, that assembled existst and may be skipping...
    if not os.path.isfile(args.assembled) or os.path.getsize(args.assembled) == 0:
        log.info("\n== Skipping processing of %s (empty file)\n" % args.assembly_type)
    else:
        log.info("\n== Processing of %s\n" % args.assembly_type)
        tmp_dir_for_corrector = os.path.join(args.output_dir, "mismatch_corrector", args.assembly_type)

        # correcting
        result_corrected_filename = os.path.join(tmp_dir_for_corrector, "corrected_contigs.fasta")

        dst_configs = os.path.join(tmp_dir_for_corrector, "configs")
        cfg_file_name = os.path.join(dst_configs, "corrector.info")

        binary_name = "spades-corrector-core"
        command = [os.path.join(args.bin_home, binary_name),
                   os.path.abspath(cfg_file_name), os.path.abspath(args.assembled)]

        log.info("\n== Running contig polishing tool: " + ' '.join(command) + "\n")
        log.info("\n== Dataset description file was created: " + cfg_file_name + "\n")
        log.info("Run: " + ' '.join(command))

        support.sys_call(command, log)

        if not os.path.isfile(result_corrected_filename):
            log.error("mismatch correction finished abnormally: %s not found!" % result_corrected_filename)

        if os.path.isfile(result_corrected_filename):
            shutil.copyfile(result_corrected_filename, args.corrected)


if __name__ == "__main__":
    main()
