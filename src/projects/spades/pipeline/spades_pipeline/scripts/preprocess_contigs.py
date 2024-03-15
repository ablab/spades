#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import logging
import sys
from os.path import abspath, dirname, realpath, join

# init python_lib_folder
python_modules_home = abspath(dirname(realpath(__file__)))
sys.path.append(join(python_modules_home, ".."))
import support


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--args_filename",
                    help="path to file with args",
                    action="store")
    parser.add_argument("--dst",
                    help="path to dst dir",
                    action="store")
    parser.add_argument("--threshold_for_breaking_additional_contigs",
                    dest="THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS",
                    type=int,
                    help="threshold for breaking additional contigs",
                    action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    # create logger
    log = logging.getLogger("Preprocess additional contigs")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    with open(args.args_filename) as f:
        lines = f.readlines()
        for gzipped, old_filename, new_filename in zip(lines[0::3], lines[1::3], lines[2::3]):
            gzipped = (gzipped.rstrip() == "True")
            old_filename = old_filename.rstrip()
            new_filename = new_filename.rstrip()

            modified, new_fasta = support.break_scaffolds(old_filename,
                                                          args.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS,
                                                          replace_char='A', gzipped=gzipped)
            log.info("== Processing additional contigs (%s): changing Ns to As and "
                     "splitting by continues (>= %d) Ns fragments (results are in %s directory)" % (
                         old_filename,
                         args.THRESHOLD_FOR_BREAKING_ADDITIONAL_CONTIGS,
                         args.dst))
            support.write_fasta(new_filename, new_fasta)


if __name__ == "__main__":
    main()
