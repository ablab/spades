#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os
import shutil
import re
import subprocess

f = open(sys.argv[1], 'r')
n50_limit = int(sys.argv[2])
mis_limit = int(sys.argv[3])

columns = map(lambda s: s.strip(), f.readline().split('\t'))
values = map(lambda s: s.strip(), f.readline().split('\t'))
n50 = int(values[columns.index("N50")])
mis = int(values[columns.index("# misassemblies")])
mismatches = float(values[columns.index("# mismatches per 100 kbp")])
indels = float(values[columns.index("# indels per 100 kbp")])
print 'N50 =', n50
print 'Misasemblies =', mis
lvl = 0
if n50 < n50_limit:
    print 'N50 is too small: less than', n50_limit
    lvl += 1
if mis > mis_limit:
    print 'Too many misassemblies: more than', mis_limit
    lvl += 2
if len(sys.argv) > 4:
    genes_limit = int(sys.argv[4])
    genes = int(values[columns.index("# genes")].split('+')[0])
    print 'full genes =', genes

    if genes < genes_limit:
	print 'Too few genes, less than', genes_limit
	lvl += 4
if len(sys.argv) > 5:
    mapped_limit = float(sys.argv[5])
    mapped = float(values[columns.index("Genome fraction (%)")])
    print 'mapped genome =', mapped, '%'

    if mapped < mapped_limit:
        print 'Too few mapped genome, less than', mapped_limit
        lvl += 8
if len(sys.argv) > 6:
    mismatches_limit = float(sys.argv[6])
    print 'mismatches per 100kb =', mismatches 

    if mismatches > mismatches_limit:
        print 'Too much mismatches, more than', mismatches_limit
        lvl += 16
if len(sys.argv) > 7:
    indels_limit = float(sys.argv[7])
    print 'indels per 100kb =', indels

    if indels > indels_limit:
        print 'Too much indels, more than', indels_limit
        lvl += 32



f.close()
sys.exit(lvl)
