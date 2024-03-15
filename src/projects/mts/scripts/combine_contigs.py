#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import os.path
from Bio import SeqIO
from common import sample_name

files = sys.argv[1:]

output = sys.stdout

for file in files:
    for seq in SeqIO.parse(file, "fasta"):
        seq.id = sample_name(file) + "-" + seq.id
        seq.description = ""
        SeqIO.write(seq, output, "fasta")
