#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import os.path
from Bio import SeqIO
from common import sample_name

replace = False

if sys.argv[1] == "-r":
    replace = True
    files = sys.argv[2:]
else:
    files = sys.argv[1:]

output = sys.stdout

for file in files:
    for seq in SeqIO.parse(file, "fasta"):
        seq_id = seq.id
        if replace:
            seq_id = seq_id.replace(",", "~")
        seq.id = sample_name(file) + "-" + seq_id
        seq.description = ""
        SeqIO.write(seq, output, "fasta")
