#!/usr/bin/env python
from __future__ import print_function

import sys
from Bio import SeqIO

if len(sys.argv) < 4:
    print("Usage:", sys.argv[0], "<contigs_file> <file with names> <output> [<operation mode>]")
    print("Operation mode is \"retain\" (default) or \"remove\"")
    sys.exit(1)

f_n = sys.argv[1]
names_f = open(sys.argv[2], "r")
names = set(l.strip() for l in names_f.readlines())
input_seq_iterator = SeqIO.parse(open(f_n, "r"), "fasta")

filtered_iterator = (record for record in input_seq_iterator \
                      if record.name in names)

if (len(sys.argv) == 5):
    if sys.argv[4] == "remove":
        filtered_iterator = (record for record in input_seq_iterator \
                      if record.name not in names)
    else:
        if sys.argv[4] != "retain":
            print("Wrong operation mode")

output_handle = open(sys.argv[3], "w")
SeqIO.write(filtered_iterator, output_handle, "fasta")
output_handle.close()
