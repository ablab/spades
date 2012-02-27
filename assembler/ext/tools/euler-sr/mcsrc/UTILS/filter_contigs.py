#!/usr/bin/env python

import sys
import os
import sequtils
import getopt
import string

True = 1
False = 0
def valid_contig(contigStr, numRepeat, numN):
    # want to make sure there aren't too many N's nor 'ACTG (non-repeat, mislabeled)'.
    nA = 0
    nC = 0
    nT = 0
    nG = 0
    nN = 0
    for i in range(0,len(contigStr)):
        if (contigStr[i] == 'A'):
            nA = nA + 1
        if (contigStr[i] == 'C'):
            nC = nC + 1
        if (contigStr[i] == 'G'):
            nG = nG + 1
        if (contigStr[i] == 'T'):
            nT = nT + 1
        if (contigStr[i] == 'N'):
            nN = nN + 1
    nACTG = nA + nC + nT + nG
    if (numRepeat > 0 and nACTG > numRepeat ):
        return False
    elif (numN > 0 and nN > numN):
        return False
    else:
        return True
    
        
        

inFileName = ""
outFileName = ""
minLength = 0
maxLength = 9999999999999L
numRepeat = -1
numN      = -1
(opts, rest) = getopt.getopt(sys.argv[1:], "i:m:M:o:r:n:")
for option,argument in opts:
    if (option == '-i'):
        inFileName = argument
    if (option == '-m'):
        minLength = int(argument)
    if (option == '-M'):
        maxLength = int(argument)
    if (option == '-o'):
        outFileName = argument
    if (option == '-r'):
        numRepeat = int(argument)
    if (option == '-n'):
        numN  = int(argument)

if (inFileName == "" or outFileName == ""):
    print "usage: filter_contigs.py -i infile -o outfile [-m min length] [-M max length ]"
    sys.exit(1)
    
inFile  = open(inFileName, 'r')
outFile = open(outFileName, 'w')

title = string.strip(inFile.readline())
count = 1

while (title != ''):
    contig, title2 = sequtils.get_contig(inFile)
    l = len(contig)
    if (l > minLength and l < maxLength):
        if (valid_contig(contig, numRepeat, numN)):
            sequtils.print_seq(outFile, contig, title + '_'+str(count)+'\n', 70)
            count = count + 1
    title = string.strip(title2)
    

