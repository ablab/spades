#!/usr/bin/env python

############################################################################
# Copyright (c) 2015-2019 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import os
import sys
from os.path import abspath, dirname, realpath, join


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result_scaffolds_filename",
                        help="path to result scaffolds filename",
                        action="store")
    parser.add_argument("--misc_dir",
                        help="path to misc dir",
                        action="store")
    parser.add_argument("--threshold_for_breaking_scaffolds",
                        dest="THRESHOLD_FOR_BREAKING_SCAFFOLDS",
                        type=int,
                        help="threshold for breaking scaffolds",
                        action="store")
    return parser.parse_args()

def main():
    args = parse_args()

    # init python_lib_folder
    python_modules_home = abspath(dirname(realpath(__file__)))
    sys.path.append(join(python_modules_home, ".."))
    import support

    if os.path.isfile(args.result_scaffolds_filename):
        if not os.path.isdir(args.misc_dir):
            os.makedirs(args.misc_dir)

        result_broken_scaffolds = os.path.join(args.misc_dir, "broken_scaffolds.fasta")
        modified, broken_scaffolds = support.break_scaffolds(args.result_scaffolds_filename,
                                                             args.THRESHOLD_FOR_BREAKING_SCAFFOLDS)
        if modified:
            support.write_fasta(result_broken_scaffolds, broken_scaffolds)


if __name__ == "__main__":
    main()
