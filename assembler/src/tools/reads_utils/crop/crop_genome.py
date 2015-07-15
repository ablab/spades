#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os

# crop genome to first n basepairs
gf = open(sys.argv[1])
n = int(sys.argv[2])
ogf = open(sys.argv[3], 'w')

for line in gf:
    ogf.write(line.strip() + ",_cropped_to_first_" + str(n) + "_nucleotides\n")
    break

cnt = 0
for line in gf:
    line = line.strip()
    if (n <= 0):
	break;
    if n < len(line):
	line = line[0:n]
    ogf.write(line)
    ogf.write('\n')
    n -= len(line)

gf.close()
ogf.close()
