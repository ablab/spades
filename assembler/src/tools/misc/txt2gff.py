#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import shutil
import sys

if len(sys.argv) != 3:
    print ("Converter for genes and operons files from .txt to .gff format")
    print ("Usage: " + sys.argv[0] + "  <input .txt file>  <keyword: either 'gene' or 'operon'>")
    sys.exit(1)

input_filename = sys.argv[1]
keyword = sys.argv[2]

if not os.path.exists(input_filename):
    print ("Specified file not found! Exiting...")
    sys.exit(1)

basename = os.path.splitext(input_filename)[0]
output_filename = basename + ".gff"

if os.path.exists(output_filename):
    print ("Output file (" + output_filename + ") already exists! Exiting...")
    sys.exit(1)

out = open(output_filename, 'w')
out.write("##gff-version	2\n")
for line in open(input_filename, 'r'):

    # EXAMPLE:
    # gi|48994873|gb|U00096.2|	1	4263805	4264884
    # gi|48994873|gb|U00096.2|	2	795085	795774 

    sections = line.split()
    if len(sections) != 4:
        continue
    ref_id  = sections[0]
    gene_id = sections[1]
    start   = int(sections[2])
    end     = int(sections[3])
    
    # EXAMPLE:    
    # ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN

    out.write(ref_id + '\t' + '.' + '\t' + keyword + '\t')
    out.write(str(min(start, end)) + '\t' + str(max(start, end)) + '\t' + '.' + '\t')
    if start < end:
        out.write('+')
    else:
        out.write('-')
    out.write('\t' + '.' + '\t' + 'ID=' + gene_id + '\n')

out.close()
print ("Conversion finished! See " + output_filename)
