#!/usr/bin/python -O

import sys
import os
import re

if len(sys.argv) < 2:
    print("Misassembled contigs getter: prints IDs of misassembled contigs or save these contigs in fasta (if input contigs file and output file specified)")
    print("Usage: " + sys.argv[0] + " <plantagota output> [input_contigs.fasta misassembled_contigs.fasta]")        
    sys.exit()

in_file = open(sys.argv[1], "r")

outFileName, ext = os.path.splitext(sys.argv[1])
outFileName += ".txt.mis"
outFile = open(outFileName, "w")

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

for contig_id in mis_contigs_ids:
    outFile.write(contig_id + '\n') 

in_file.close()
outFile.close()

