#!/usr/bin/python

#add some  contamination for meta-genomic simulation

import sys
import os
import random

if len(sys.argv) != 5:
	print("Usage: " + sys.argv[0] + " <part> <contamination_reads> <main_reads> <output_file>")	
	print ("add 1/<part> reads from contamination_reads to main reads")
	sys.exit()

contFileName = sys.argv[2]
contFile = open(contFileName, "r")
mainFileName = sys.argv[3]
mainFile = open(mainFileName, "r")
outFileName = sys.argv[4]
outFile = open(outFileName, "w")

threshold = 1/float(sys.argv[1])



mainName, mainext = os.path.splitext(mainFileName)
contName, context = os.path.splitext(contFileName)

#outFile = open(mainName + "+"+ str(threshold) + "* contamunation"  + ext, "w") 
line = contFile.readline()
while 1:
	if not line:
		break
	out = random.random() < threshold
	for i in range(8):
		if out:
			outFile.write(line)
		line = contFile.readline()



line = mainFile.readline()

while line:
	outFile.write(line)
	line = mainFile.readline()


#c_len, contig = read_contig(inFile)
#while contig is not None:
#	if c_len <= threshold:
#		fileS.write(contig)
#	else:
#		fileL.write(contig)
#		
#	c_len, contig = read_contig(inFile)


contFile.close()
mainFile.close()
outFile.close()


