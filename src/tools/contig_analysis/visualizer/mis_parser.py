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

basename, ext = os.path.splitext(sys.argv[1])
out_IDs_filename = basename + ".txt.mis"
out_aligns_filename = basename + ".aligns"
out_IDs_file = open(out_IDs_filename, "w")
out_aligns_file = open(out_aligns_filename, "w")

mis_contigs_ids = []
mis_contigs_aligns = []

#skipping prologue
for line in in_file:
    if line.startswith("Analyzing contigs..."):
        break

# main part of plantagora output
cur_contig_id = ""
cur_contig_align = ""
was_extensive_mis = False
for line in in_file:
    if line.startswith("Analyzing coverage..."):
        break

    if line.startswith("	CONTIG:"):
        cur_contig_id = line.split("	CONTIG:")[1].strip()
        if was_extensive_mis:
            was_extensive_mis = False
            mis_contigs_aligns.append(cur_contig_align)
        cur_contig_align = ""
    cur_contig_align += line

    if (line.find("Extensive misassembly") != -1) and (cur_contig_id != ""):
        was_extensive_mis = True
        mis_contigs_ids.append(cur_contig_id.split()[0])
        cur_contig_id = ""

if was_extensive_mis:
    mis_contigs_aligns.append(cur_contig_align)
            
# printing IDs of misassembled contigs

print "IDs of misassembled contigs saved in " + out_IDs_filename
for contig_id in mis_contigs_ids:
    out_IDs_file.write(contig_id + '\n') 

print "Alignments of misassembled contigs saved in " + out_aligns_filename
for contig_align in mis_contigs_aligns:
    out_aligns_file.write(contig_align) 

in_file.close()
out_IDs_file.close()
out_aligns_file.close()
