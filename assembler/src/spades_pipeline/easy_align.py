#!/usr/bin/env python
import os
import sys

pipeline_modules_home = 'src/spades_pipeline/'#os.path.dirname(os.path.realpath(__file__)
sys.path.append(os.path.join(pipeline_modules_home, "common"))
sys.path.append(os.path.join(pipeline_modules_home, "truspades"))

import logging
import alignment

#import spades_init
#spades_init.init()


def align(contigs_file, left, right, out_dir, threads):
    #logging
    log = logging.getLogger('reference_construction')
    log.setLevel(logging.INFO)
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.INFO)
    log.addHandler(console)
    #logging

    bwa_command='bin/spades-bwa'
    index = alignment.index_bwa(bwa_command, log, contigs_file, os.path.join(out_dir, "bwa_index"), "bwtsw")
    index = os.path.join(out_dir, "bwa_index", "index")
    sam = alignment.align_bwa_pe_lib(bwa_command, index, left, right, os.path.join(out_dir, "align"), log, threads)
    #index_bwa(command, log, reference, work_dir, algorithm = "is"):
    #align_bwa_pe_lib(command, work_dir + "/index", reads_file1, reads_file2, work_dir, log, threads = 1):

if __name__ == '__main__':

    if len(sys.argv) < 5:
        sys.stderr.write("Usage: %s <contigs> <left_reads> <right_reads> <out_dir> [threads = 8]\n" % sys.argv[0])
        exit(1)

    threads = 8
    if len(sys.argv) >= 6:
        threads = int(sys.argv[5])

    align(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], threads)
