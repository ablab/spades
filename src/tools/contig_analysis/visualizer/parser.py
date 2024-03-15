#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os

if len(sys.argv) != 2:
	print("Usage: " + sys.argv[0] + " <plantagota output>")	
	sys.exit()

inFileName = sys.argv[1]

outFileName, ext = os.path.splitext(inFileName)
outFileName += ".txt"
outFile = open(outFileName, "w")

misFileName, ext = os.path.splitext(inFileName)
misFileName += ".txt.mis"
misFile = open(misFileName, "w")


inFile = open(inFileName)
mis_contigs_ids = []

#skipping prologue
for line in inFile:
	if line.startswith("Analyzing contigs..."):
		break

# main part of plantagora output
cur_contig_id = ""
for line in inFile:
	if line.startswith("	CONTIG:"):
		cur_contig_id = line.split("	CONTIG:")[1].strip()
	if (line.find("Extensive misassembly") != -1) and (cur_contig_id != ""):
		mis_contigs_ids.append(cur_contig_id.split()[0])
		cur_contig_id = ""
	if line.startswith("Analyzing coverage..."):
		break

mis_contigs = {}

inFile = open(inFileName)
for line in inFile:
	parse = line.strip().split(' ')
	if parse[0] == "Align":
		forward = " +"
		if (int(parse[2]) - int(parse[3]))*(int(parse[5]) - int(parse[6])) < 0:
			forward = " -"
		cid = parse[4]
		outFile.write(parse[2] + " " + parse[3] + " " + cid + forward + " " + str(min(int(parse[5]),int(parse[6]))) + "\n")

		if cid in mis_contigs_ids:
			if cid not in mis_contigs:
				mis_contigs[cid] = ""
			mis_contigs[cid] += cid + " " + parse[2] + " " + parse[3] + " " + forward + " " + str(min(int(parse[5]),int(parse[6]))) + "\n"

for key in mis_contigs:
	misFile.write(mis_contigs[key])

print(outFileName)

outFile.close()
inFile.close()
misFile.close()
		
	
