import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio import SeqIO
import sys
import os

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
len_seq = sorted([len(record.seq) for record in input_seq_iterator])

buckets = {2000: 0, 5000: 0, 10000: 0}

for l in len_seq:
	if l < 5000:
		buckets[2000] += 1
	elif l < 10000:
		buckets[5000] += 1
	else:
		buckets[10000] += 1

print buckets

