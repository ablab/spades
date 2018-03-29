#!/usr/bin/python

#add some  contamination for meta-genomic simulation

import sys
import os
import random
from os import listdir
from os.path import isfile, join
import fastaparser
if len(sys.argv) != 3:
    print "Usage: " + sys.argv[0] + " <directory with fastas> <kmer> "
    exit()
dir = sys.argv[1]
kmer = int(sys.argv[2])
#to avoid problems with repeats caused by circular structure
back = 400
if (kmer < back):
    back -= kmer
else:
    back = 1
unique = 0
non_unique = 0
for file in os.listdir(dir):
    arr = file.split('.')
    if arr[-1] != "fasta":
        continue
    out_dir = join(dir, arr[0])
    contigs = fastaparser.read_fasta(join(dir,file))
    tr_file = join(dir, (arr[0] + "_tr.fsa"));
    for contig in contigs:
        contig_n = [contig[0], contig[1][:-back]]
        os.system("rm " + tr_file)
        fastaparser.write_fasta_to_file(tr_file, [contig_n])
    jf_file = join(dir,arr[0]+ ".jf")
    jf_stats = join(dir, arr[0] + ".stats")
    os.system("jellyfish count " + tr_file +  " -m " + str(kmer)+" -C -c 2 -s 3000000 -o "+ jf_file)
    os.system("jellyfish stats " + jf_file + " > " + jf_stats)
    with open (jf_stats,"r") as stats:
        lines = stats.readlines()
        max_count = int(lines[3].split()[1])
        length = int(lines[2].split()[1])
#        print " max multiplicity " + str(max_count)
        kmer_repeat = ""
        if max_count == 1:
            kmer_repeat = "NO"
            unique += 1
        else:
            kmer_repeat = "YES"
            non_unique += 1

        print os.path.basename(jf_stats)+ " \t " + str(length + back + kmer - 1) + " \t " +  kmer_repeat
#        print os.path.basename(jf_stats)+ " len: " + str(length + back + kmer - 1) + " max rep multiplicity: " +  str(max_count)
print "With repeat size " + str(kmer) + " no repeats: " +  str(unique) + " repeats:  " + str(non_unique)


