#!/usr/bin/env python
from __future__ import print_function
try:
    from itertools import izip as zip
except ImportError:
    pass

import re
import argparse
import os
import sys
from shutil import copyfile

parser = argparse.ArgumentParser(description="Binner input formatter")
parser.add_argument("--prefix", "-o", type=str, help="Prefix path to save files for MaxBin")
parser.add_argument("--profile", "-p", type=str, help="Path to concoct profile")

args = parser.parse_args()

#copyfile(args.profile, args.file)
out_list = open(args.prefix + "_list.txt", "w")
prof = open(args.profile, "r")
samples_files = []
first = True
for ln in prof.readlines():
    if first:
        first = False
        lst = ln.strip().split("\t")
        for ind in range(1, len(lst)):
            samples_files.append(open(args.prefix + "_" + str(ind), "w"))
            print(args.prefix + "_" + str(ind), file=out_list)
    lst = ln.strip().split("\t")
    contig_name = lst[0]
    for ind in range(len(samples_files)):
        samples_files[ind].write("\t".join([contig_name, lst[ind + 1] ]) + "\n")

for sf in samples_files:
    sf.close()

out_list.close()
prof.close()
