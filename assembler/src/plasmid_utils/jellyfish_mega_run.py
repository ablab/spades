#!/usr/bin/python

#add some  contamination for meta-genomic simulation

import sys
import os
import random
from os import listdir
from os.path import isfile, join
import fastaparser
dir = sys.argv[1]
unique = 0
non_unique = 0
for file in os.listdir(dir):
    arr = file.split('.')
    if arr[-1] != "fasta":
        continue
    out_dir = join(dir, arr[0])
    contigs = fastaparser.read_fasta(join(dir,file))
    tr_file = join(dir, (arr[0] + "_tr.fasta"));
    for contig in contigs:
        contig_n = [contig[0], contig[1][:-345]]
        fastaparser.write_fasta_to_file(tr_file, [contig_n])
    jf_file = join(dir,arr[0]+ ".jf")
    jf_stats = join(dir, arr[0] + ".stats")
    os.system("jellyfish count " + tr_file +  " -m 55 -C -c 2 -s 3000000 -o "+ jf_file)
    os.system("jellyfish stats " + jf_file + " > " + jf_stats)
    with open (jf_stats,"r") as stats:
        lines = stats.readlines()
        max_count = int(lines[3].split()[1])
        print max_count
        if max_count == 1:
            unique += 1
        else:
            non_unique += 1
        print str(unique) + " " + str(non_unique)

