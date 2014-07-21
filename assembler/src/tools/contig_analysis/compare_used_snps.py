#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys

def read_snps(filename):
    """
        Returns set of genome snps positions
    """
    
    res_snps = dict()
    for line in open(filename):
        arr = line.split()
        pos = int(arr[2])
        res_snps[pos] = arr[1] + ' ' + arr[2] + ' ' + arr[5]
    return res_snps





if len(sys.argv) != 3:
	print('Used snps comparator usage: ' + sys.argv[0] + ' <file1> <file2>')
	sys.exit(0)



first_snps = read_snps(sys.argv[1])
second_snps = read_snps(sys.argv[2])
print("SNPs present in " + sys.argv[1] + " not present in " + sys.argv[2] + ':')
for snp in first_snps:
    if not snp in second_snps:
        print  first_snps[snp]

print("\n SNPs present in " + sys.argv[2] + " not present in " + sys.argv[1] + ':')
for snp in second_snps:
    if not snp in first_snps:
        print  second_snps[snp]