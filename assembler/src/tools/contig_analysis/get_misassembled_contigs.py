#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import re

if len(sys.argv) < 2:
    print("Misassembled contigs getter: prints IDs of misassembled contigs or save these contigs in fasta (if input contigs file and output file specified)")
    print("Usage: " + sys.argv[0] + " <plantagota output> [input_contigs.fasta misassembled_contigs.fasta]")        
    sys.exit()

in_file = open(sys.argv[1], "r")

mis_contigs_ids = []

#skipping prologue
for line in in_file:
    if line.startswith("Analyzing contigs..."):
        break

# main part of plantagora output
cur_contig_id = ""
for line in in_file:
    if line.startswith("	CONTIG:"):
        cur_contig_id = line.split("	CONTIG:")[1].strip()
    if (line.find("Extensive misassembly") != -1) and (cur_contig_id != ""):
        mis_contigs_ids.append(cur_contig_id.split()[0])
        cur_contig_id = ""
    if line.startswith("Analyzing coverage..."):
        break
            
# printing IDs of misassembled contigs
print("Misassembled contigs:")
for contig_id in mis_contigs_ids:
    print(contig_id) 

in_file.close()

if (len(sys.argv) == 4):
    import fastaparser
    input_contigs = fastaparser.read_fasta(sys.argv[2])
    mis_contigs = open(sys.argv[3], "w")    

    for (name, seq) in input_contigs:
        corr_name = re.sub(r'\W', '', re.sub(r'\s', '_', name))

        if mis_contigs_ids.count(corr_name) != 0:
            mis_contigs.write(name + '\n')
            for i in xrange(0, len(seq), 60):
                mis_contigs.write(seq[i:i+60] + '\n')

    mis_contigs.close()
