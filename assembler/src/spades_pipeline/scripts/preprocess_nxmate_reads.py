#!/usr/bin/env python

############################################################################
# Copyright (c) 2019 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import argparse
import logging
import sys
from os.path import abspath, dirname, realpath, join, isfile

# init python_lib_folder
python_modules_home = abspath(dirname(realpath(__file__)))
sys.path.append(join(python_modules_home, ".."))

if isfile(join(python_modules_home, "../../../spades_init.py")):
    sys.path.append(join(python_modules_home, "../../.."))
else:
    sys.path.append(join(python_modules_home, "../../../../bin/"))

import lucigen_nxmate, support


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--args_filename",
                        help="path to file with args",
                        action="store")
    parser.add_argument("--dst",
                        help="path to dst dir",
                        action="store")
    parser.add_argument("--threads",
                        help="number of threads",
                        default=16,
                        type=int,
                        action="store")
    return parser.parse_args()

def main():
    args = parse_args()

    # create logger
    log = logging.getLogger("Preprocess Lucigen NxMate reads")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    try:
        with open(args.args_filename) as f:
            lines = f.readlines()
            for infile1, infile2 in zip(lines[0::2], lines[1::2]):
                lucigen_nxmate.process_reads(infile1, infile2, args.dst, log, args.threads)
    except ImportError:
        support.error("can't process Lucigen NxMate reads! lucigen_nxmate.py is missing!", log)


if __name__ == "__main__":
    main()
