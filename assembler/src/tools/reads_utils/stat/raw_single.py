#!/usr/bin/python -O

#Create raw file for single reads using chromosome map file (optional)
#Chromosome map file should contain choromosome name (as used in bowtie) and its start position in genome, one chromosome per line

import sys

def raw_single(in_filename, out_filename, chrsm_file = None):
    inFile = open(in_filename)
    outFile = open(out_filename, 'w')

    wChrs = (chrsm_file != None)

    chrs = {" ":0}
    if wChrs:	 
        chrsF = open(chrsm_file)
        for line in chrsF:
            chrs[line.split()[0]] =	int(line.split()[1])
        chrsF.close()

    delim = '/'

    for line in inFile:
        chr1 = line.split('\t',5)[2]
        pos1 = int(line.split('\t', 5)[3])
        len1 = len(line.split('\t', 5)[4])

        addPos = 0
        if wChrs:
            addPos = chrs[chr1] 	
		
        outFile.write(str(pos1 + addPos) + ' ' + str(len1) + '\n')

    inFile.close()
    outFile.close()

def main():

    if len(sys.argv) != 4 and len(sys.argv) != 3:
        print("Usage: <bowtie log file> <output> [chromosomes file]");
        exit(0)

    if (len(sys.argv) == 4):
        raw_single(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        raw_single(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
