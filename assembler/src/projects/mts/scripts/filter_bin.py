#!/usr/bin/env python
from __future__ import print_function

import sys
from Bio import SeqIO
import common

def print_usage():
        print("Usage: filter_bins.py <contigs> <binning info> <bin name>")

contigs = sys.argv[1]
binning = common.load_annotation(sys.argv[2], False)
bin_name = sys.argv[3]

for seq in SeqIO.parse(contigs, "fasta"):
    if bin_name in binning.get(seq.id, set()):
        SeqIO.write(seq, sys.stdout, "fasta")
