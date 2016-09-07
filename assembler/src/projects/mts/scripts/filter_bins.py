#!/usr/bin/python
from __future__ import print_function

import sys
import common

def print_usage():
    print("Usage: filter_bins.py <annotation> <contigs>")

binning = common.load_annotation(sys.argv[1])

bins = set(bin for contig in open(sys.argv[2]) for bin in binning.get(contig, []))
for bin in bins:
    print(bin)
