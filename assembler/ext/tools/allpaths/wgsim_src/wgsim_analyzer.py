#!/usr/bin/python

import sys
import os
import shutil
import re

#################
#
# 
# @NODE_1_length_42021_cov_807.861328_26499_26676_4:0:0_3:0:0_0/1
# GGCGGGCAATGATTGGGTTCCTACTTCGATTACCGCTTATCTGGCGGGCGGGATGTTTTTACAATGGCTGCTGGGGCCGCTGTCGGATCGTATTGGTCGC
# +
# 2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
# @NODE_1_length_42021_cov_807.861328_23575_23760_3:0:0_2:0:0_1/1
# GATGATGCACGGTTTTTATAGCCTGGGCACGCTGGCCGGCGCAGGTGTCGGGATGGCACTGACGGCCTTTGGCGTTCCGGCCACGGTGCACATTTTATTT
# +
# 2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
# 
# 
#################

def read_block(in_file):
    block = []
    for i in range(4):
        line = in_file.readline()
        if not line:
            return block
        block.append(line)
    return block   

def reverse(seq_str):
    compl = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    nucls = set(compl.keys())

    seq_list = list(seq_str)
    for i in range(0, len(seq_list)):
        if seq_list[i] in nucls:
            seq_list[i] = compl[seq_list[i]]

    return "".join(seq_list)

########################################################################	

if len(sys.argv) < 3:
	print 'Script for reverse complement fq reads'
	print 'Usage:', sys.argv[0], 'reads1.fq reads2.fq out.mf'
	sys.exit(0)

reads1 = sys.argv[1]
reads2 = sys.argv[2]
out    = sys.argv[3]

reads1_file = open(reads1, 'r')
reads2_file = open(reads2, 'r')
out_file    = open(out, 'w')

block1 = read_block(reads1_file)
block2 = read_block(reads2_file)

while (len(block1) == 4) and (len(block2) == 4):
    if block1[0].split('/')[0] == block2[0].split('/')[0]:
        out_file.write(block1[0].split('/')[0] + '/1 and /2\n')
    out_file.write(block1[1])
#    out_file.write(block2[1])
    out_file.write(reverse(block2[1]))
    
    block1 = read_block(reads1_file)
    block2 = read_block(reads2_file)
     
reads1_file.close()
reads2_file.close()
out_file.close()
