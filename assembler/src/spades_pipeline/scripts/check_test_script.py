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
import sys
from os.path import abspath, dirname, realpath, join


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode",
                        choices=["common", "truseq", "rna", "plasmid"],
                        help="running mode",
                        action="store")
    parser.add_argument("--truseq_long_reads_file",
                        help="path to truseq long reads",
                        action="store")
    parser.add_argument("--result_transcripts_filename",
                        help="path to file with result transcripts",
                        action="store")
    parser.add_argument("--result_contigs_filename",
                        help="path to file with result contigs",
                        action="store")
    parser.add_argument("--result_scaffolds_filename",
                        help="path to file with result scaffolds",
                        action="store")
    return parser.parse_args()


def main():
    args = parse_args()

    # init python_lib_folder
    python_modules_home = abspath(dirname(realpath(__file__)))
    sys.path.append(join(python_modules_home, ".."))
    import support

    # create logger
    log = logging.getLogger("Check test")
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(message)s"))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)

    if args.mode == "truseq":
        if not os.path.isfile(args.truseq_long_reads_file):
            support.error("TEST FAILED: %s does not exist!" % args.truseq_long_reads_file)
    elif args.mode == "rna":
        if not os.path.isfile(args.result_transcripts_filename):
            support.error("TEST FAILED: %s does not exist!" % args.result_transcripts_filename)
    else:
        for result_filename in [args.result_contigs_filename, args.result_scaffolds_filename]:
            if os.path.isfile(result_filename):
                result_fasta = list(support.read_fasta(result_filename))
                # correctness check: should be one contig of length 1000 bp
                correct_number = 1
                if args.mode == "plasmid" or args.mode == "metaplasmid":
                    correct_length = 9689
                else:
                    correct_length = 1000
                if not len(result_fasta):
                    support.error("TEST FAILED: %s does not contain contigs!" % result_filename)
                elif len(result_fasta) > correct_number:
                    support.error("TEST FAILED: %s contains more than %d contig (%d)!" %
                                  (result_filename, correct_number, len(result_fasta)))
                elif len(result_fasta[0][1]) != correct_length:
                    if len(result_fasta[0][1]) > correct_length:
                        relation = "more"
                    else:
                        relation = "less"
                    support.error("TEST FAILED: %s contains %s than %d bp (%d bp)!" %
                                  (result_filename, relation, correct_length, len(result_fasta[0][1])))
            else:
                support.error("TEST FAILED: %s does not exist!" % result_filename)
    log.info("\n========= TEST PASSED CORRECTLY.")


if __name__ == "__main__":
    main()
