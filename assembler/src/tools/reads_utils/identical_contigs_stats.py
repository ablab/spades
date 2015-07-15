#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate statistics for reads by reference file

import argparse
import sys
import os

if len(sys.argv) != 2:
        print("Usage: " + sys.argv[0] + "  <contigs_file>")     
	print("identical contigs names outputs to repeated_contigs_names.txt file")
        sys.exit()

contFileName = sys.argv[1]
contFile = open(contFileName, "r")
outFileName = "repeated_contigs_names.txt";
outFile = open(outFileName, "w")


line = contFile.readline()
contigs = set([])
exact_repeats = 0;
total_contigs = 0;
while 1:
	if not line:
		break
#	print(line)
	if (line[0] != '>'):
		print("something wrong" + str(line[1]))
		break
	contig = ""
	contig_name = line;
	line = contFile.readline()
#	print(line);
	while ((line[0] !='>')):
		contig += line;
		line = contFile.readline()		
	        if not line:
        	        break

	if (contig in contigs):
		outFile.write(contig_name)
		exact_repeats += 1
	total_contigs += 1
	contigs.add(contig);	
print ("Exact repeats: " + str(exact_repeats) + "; total contigs: " + str(total_contigs));	

