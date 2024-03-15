#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os

DATA_PATH = '../'

max = 100000
if len(sys.argv) >= 2:
    max = int(sys.argv[1])

# crop reference genome to first 'max' basepairs
filename1 = DATA_PATH + '/ref.fasta' # input genome
filename2 = 'ref.' + str(max) + '.fasta' # output genome (cropped)
f1 = open(filename1)
f2 = open(filename2, 'w')
header = f1.readline()
f2.write(header)
cnt = 0 # because of the '>' in first line
for line in f1:
    line = line.strip()
    if max - cnt < len(line):
        line = line[0:max-cnt]
    f2.write(line + '\n')
    cnt += len(line)
    if cnt >= max:
        break
f1.close()
f2.close()

# build bowtie index
index_name = 'bowtie_index.'+str(max)
os.system('bowtie-build ' + filename2 + ' ' + index_name)

# align reads using bowtie
os.system('bowtie ' + index_name + ' -1 ' + DATA_PATH + '/std_left.cor.fastq -2 ' + DATA_PATH + '/std_right.cor.fastq --al std.first' + str(max) + '.fastq')
# use -X for maxins length

# gzip resulting files
# os.system('gzip frag.first' + str(max) + '_1.fastq')
# os.system('gzip frag.first' + str(max) + '_2.fastq')
