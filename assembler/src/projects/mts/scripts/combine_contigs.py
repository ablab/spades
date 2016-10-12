#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import re
from Bio import SeqIO

replace = False

if sys.argv[1] == "-r":
    replace = True
    files = sys.argv[2:]
else:
    files = sys.argv[1:]

sample_re = re.compile("sample\d+")

output = sys.stdout

for file in files:
    sample = sample_re.search(file).group(0)
    for seq in SeqIO.parse(file, "fasta"):
        seq_id = seq.id
        if replace:
            seq_id = seq_id.replace(",", "~")
        seq.id = sample + "-" + seq_id
        seq.description = ""
        SeqIO.write(seq, output, "fasta")
