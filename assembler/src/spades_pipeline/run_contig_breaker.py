#!/usr/bin/env python

import sys
import break_by_coverage
from Bio import SeqIO

if len(sys.argv) < 4:
    sys.stderr.write("Usage: %s <contigs> <sam_file> <output_filename>\n" % sys.argv[0])
    exit(1)

contigs_file = sys.argv[1]
sam_file = sys.argv[2]
output_file = sys.argv[3]

coverage_breaker = break_by_coverage.ContigBreaker(contigs_file, sam_file, 100, 50)
contigs = list(SeqIO.parse(open(contigs_file, "rU"), "fasta"))
output = open(output_file, "w")
for contig in contigs:
    for subcontig in coverage_breaker.Break(contig):
        SeqIO.write(subsubsubcontig, output, "fasta")
output.close()
