#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate coverage from raw file

import sys


def main():

    if len(sys.argv) < 4:
	    print("Usage: <coverage file> <output> <bin size> ");
	    exit(0)

    inFile = open(sys.argv[1])
    outFile = open(sys.argv[2], 'w')

    bin = int(sys.argv[3])
    
    hist = []
    for line in inFile:
            hist.append(int(line.split()[1]))
	    
    maxLen = len(hist)
    print(maxLen)
    newHist = [0 for i in range((maxLen + 1) / bin + 2)]
    zeroList = []

    for i in range(maxLen):
            hpos = int(i/bin)
	    newHist[hpos] += hist[i]

            if hist[i] == 0:
                    zeroList.append(hpos)

    zeroSet = set(zeroList)

    for i in range((maxLen + 1) / bin + 1):
	    outFile.write(str(i) + ' ' + str(newHist[i] / bin) + '\n')
            if i in zeroSet:
                    outFile.write(str(i) + ' 0\n')

    inFile.close()
    outFile.close()


if __name__ == '__main__':
    main()
