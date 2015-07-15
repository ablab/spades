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
import sets

def read_contigs(filename):
    res_seq = {}

    seq =''
    cont_name = ''
    for line in open(filename):
        if line[0] == '>':
            if cont_name != '':
                res_seq[cont_name] = seq
            seq = ''
            cont_name = line[1:33].strip()
        else:
            seq += line.strip()
    if cont_name != '':
        res_seq[cont_name] = seq
    return res_seq


def write_fasta(data, filename):
    outFile = open(filename, 'w')
    for seq in data:
        outFile.write('>' + seq[0].strip() + '\n' );
        i = 0
        while i < len(seq[1]):
            outFile.write(seq[1][i:i+60] + '\n')
            i += 60
    outFile.close()



if len(sys.argv) < 1  :
	print("Usage: " + sys.argv[0] + " <first_snps_file> <second_snps_file>")
	sys.exit()

infile = open(sys.argv[1],'r');
outfile = open(sys.argv[1]+"_processed",'w');

first = 1
for line in infile:

    arr = (line.strip()).split('\t')
    if len(arr) > 4:

        new_line = ""
        interesting_indexes = [0,11,15,19,21,22,23,24]
        for i in interesting_indexes:
            tmp = arr[i].split('+')[0]
            if i == 15 and not first:
                tmp = int(tmp) + int(arr[i+1].split('+')[0])
	    new_line += str(tmp) 
	    if i != 24:
		new_line += ' &\t'
            else:
                new_line += '\\\\\n'
        first = 0
        outfile.write(new_line)
    else:
        outfile.write(line)
