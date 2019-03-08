#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys

pipeline_modules_home = 'src/spades_pipeline/'  # os.path.dirname(os.path.realpath(__file__)
sys.path.append(os.path.join(pipeline_modules_home, "common"))
sys.path.append(os.path.join(pipeline_modules_home, "truspades"))

# import alignment
import sam_parser
import break_by_coverage
import SeqIO


def break_contigs(contigs_file, sam_file, output_file):
    contigs = list(SeqIO.parse(open(contigs_file, "rU"), "fasta"))
    # sam = sam_parser.SamChain([sam_parser.Samfile(sam_file) for sam_file in sam_files])
    sam = sam_parser.Samfile(sam_file)
    # last two arguments: K, min0 stretch length to break
    coverage_breaker = break_by_coverage.ContigBreaker(contigs, sam, 100, 50)
    coverage_breaker.OutputBroken(output_file)
    # contigs = list(SeqIO.parse(open(contigs_file, "rU"), "fasta"))
    # output = open(output_file, "w")
    # for contig in contigs:
    #    for subcontig in coverage_breaker.Break(contig):
    #        SeqIO.write(subcontig, output, "fasta")
    # output.close()


if __name__ == '__main__':

    if len(sys.argv) < 4:
        sys.stderr.write("Usage: %s <contigs> <sam_file> <output_filename>\n" % sys.argv[0])
        exit(1)

    contigs_file = sys.argv[1]
    sam_file = sys.argv[2]
    output_file = sys.argv[3]
    break_contigs(contigs_file, sam_file, output_file);
