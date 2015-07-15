#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#add some  contamination for meta-genomic simulation

import sys
import os
import random

if len(sys.argv) < 4 or (len(sys.argv) %2 != 0) :
	print("Usage: " + sys.argv[0] + " {<part> <contamination_reads>}   <output_file>")
	print ("add 1/<part> reads from contamination_reads to main reads")
	sys.exit()
len_arg = len(sys.argv)
file_names = []
files = []
parts = []
ds = int(len(sys.argv)/2 -1)
print len(sys.argv)
for i in range(0, ds ):
    print i;
    file_names.append(sys.argv[i * 2 + 2])
    parts.append(sys.argv[i * 2 + 1])
    files.append(open(file_names[i]))

outFileName = sys.argv[len_arg-1];
outFile = open(outFileName, "w")

threshold = 1/float(sys.argv[1])



#outFile = open(mainName + "+"+ str(threshold) + "* contamunation"  + ext, "w") 
for i in range(0, ds):
    print " mixing "+str(i) +"   of "  + str(len(files))
    line = (files[i]).readline()
    while 1:
        if not line:
            break
        out = random.random() < (1/float(parts[i]))
        for j in range(8):
            if out:
                outFile.write(line)
            line = (files[i]).readline()

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


