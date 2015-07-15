#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Create raw file for single reads using chromosome map file (optional)
#Chromosome map file should contain choromosome name (as used in bowtie) and its start position in genome, one chromosome per line

import sys

def raw_single(in_filename, out_filename):
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')

    delim = '/'

    for line in inFile:
        if line[0] == '@':
            continue

        linesplt = line.split('\t')
        pos1 = int(linesplt[3])
        len1 = len(linesplt[9].strip())
        flag = int(linesplt[1])

        if not (flag & 4): # The read has no reported alignments
            outFile.write(str(pos1) + ' ' + str(len1) + '\n')

    inFile.close()
    outFile.close()

def main():

    if len(sys.argv) != 3:
        print("Usage: <bowtie log file> <output>");
        exit(0)

    raw_single(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
