#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Make three (fr,rf,ff) raw files from one bowtie log file using chromosome map file (optional)
#Chromosome map file should contain choromosome name (as used in bowtie) and its start position in genome, one chromosome per line

import sys
import os

if len(sys.argv) != 3 and len(sys.argv) != 2:
	print("Usage: " + sys.argv[0] + " <bowtie log> [chromosomes file]")	
	sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName, 'r')

wChrs = (len(sys.argv) == 3)

chrs = {" ":0}
if wChrs:	 
	chrsF = open(sys.argv[3])
	for line in chrsF:
		chrs[line.split()[0]] =	int(line.split()[1])
	chrsF.close()

fName, ext = os.path.splitext(inFileName)
frFile = open(fName + "_FR" + ext, "w")
rfFile = open(fName + "_RF" + ext, "w")
ffFile = open(fName + "_FF" + ext, "w")
 
delim = '/'

ffC = 0
frC = 0
rfC = 0
aaC = 0

prevLine = ""
for line in inFile:
	if (prevLine == ""):
		prevLine = line
		continue

	if (line.split(delim, 1)[0] == prevLine.split(delim, 1)[0]):
		l1 = line.split('\t', 5)
		l2 = prevLine.split('\t', 5)

		or2 = l1[1]
		or1 = l2[1]
		chr2 = l1[2]
		chr1 = l2[2]
		pos2 = int(l1[3])
		pos1 = int(l2[3])
               	len2 = len(l1[4].strip())
               	len1 = len(l2[4].strip())

		if chr1 == chr2:
			addPos = 0
			if wChrs:
				addPos = chrs[chr1] 

			if (or1 == or2):
				ffC += 1
				if (pos1 < pos2):
					ffFile.write(str(pos1 + addPos) + ' ' + str(len1) + '\n' + str(pos2 + addPos) + ' ' + str(len2) + '\n')
				else:
					ffFile.write(str(pos2 + addPos) + ' ' + str(len2) + '\n' + str(pos1 + addPos) + ' ' + str(len1) + '\n')

			elif (pos1 <= pos2 and or1 == '+'):
				frC += 1
				frFile.write(str(pos1 + addPos) + ' ' + str(len1) + '\n' + str(pos2 + addPos) + ' ' + str(len2) + '\n')

			elif (pos2 < pos1 and or2 == '+'):
		                frC += 1
		                frFile.write(str(pos2 + addPos) + ' ' + str(len2) + '\n' + str(pos1 + addPos) + ' ' + str(len1) + '\n')

		        elif (pos1 <= pos2 and or1 == '-'):
		                rfC += 1
		                rfFile.write(str(pos1 + addPos) + ' ' + str(len1) + '\n' + str(pos2 + addPos) + ' ' + str(len2) + '\n')

		        elif (pos2 < pos1 and or2 == '-'):
		                rfC += 1
		                rfFile.write(str(pos2 + addPos) + ' ' + str(len2) + '\n' + str(pos1 + addPos) + ' ' + str(len1) + '\n')

			else:
				print('Something wrong: '+str(pos1)+' '+str(len1)+' '+or1+' '+str(pos2)+' '+str(len2)+' '+or2+'\n')
		else:
			aaC += 1

	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])
		exit(0)

	prevLine = ""


print(str(frC)+' '+str(rfC)+' '+str(ffC)+' '+str(aaC)+'\n')

inFile.close()
frFile.close()
rfFile.close()
ffFile.close()
