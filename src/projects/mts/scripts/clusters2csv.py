#!/usr/bin/env python
from __future__ import print_function
import sys

from Bio import SeqIO

from os import listdir
from os.path import isfile, join

if len(sys.argv) < 3:
    print("Usage: %s <cluster directory> <output> " % sys.argv[0])
    sys.exit(1)

path = sys.argv[1]

with open(sys.argv[2], "w") as output:
    for f in listdir(path):
        if isfile(join(path, f)) and f.endswith("fna"):
            cluster = f.split(".")[0].split("_")[-1]
            record_dict = SeqIO.to_dict(SeqIO.parse(join(path, f), "fasta"))
            for k in record_dict.keys():
                print(str(k) + "," + str(cluster), file=output)
