#!/usr/bin/python

#add some  contamination for meta-genomic simulation

import sys
import os
import random

if len(sys.argv) < 6 or (len(sys.argv) %2 != 0) :
	print("Usage: " + sys.argv[0] + " {<part> <contamination_reads>}   <output_file>")
	print ("add 1/<part> reads from contamination_reads to main reads")
	sys.exit()
len_arg = len(sys.argv)
file_names = []
files = []
parts = []
for i in range(0, len/2 - 1 ):
    file_names.append(sys.argv[i + 2])
    parts.append(sys.argv[i + 1])
    files.append(open(file_names[i + 2]))

outFileName = sys.argv[len_arg-1];
outFile = open(outFileName, "w")

threshold = 1/float(sys.argv[1])



mainName, mainext = os.path.splitext(mainFileName)
contName, context = os.path.splitext(contFileName)

#outFile = open(mainName + "+"+ str(threshold) + "* contamunation"  + ext, "w") 
line = contFile.readline()
for i in range(len(file_names)):

    while 1:
        if not line:
            break
        out = random.random() < (1/float(parts[i]))
        for i in range(8):
            if out:
                outFile.write(line)
            line = files.readline()

    files[i].close()



#c_len, contig = read_contig(inFile)
#while contig is not None:
#	if c_len <= threshold:
#		fileS.write(contig)
#	else:
#		fileL.write(contig)
#		
#	c_len, contig = read_contig(inFile)


outFile.close()


