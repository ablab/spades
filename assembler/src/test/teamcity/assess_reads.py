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
aligned_limit = float(sys.argv[2])
mapped_limit = float(sys.argv[3])

columns = map(lambda s: s.strip(), f.readline().split('\t'))
values = map(lambda s: s.strip(), f.readline().split('\t'))
aligned = values[columns.index("Aligned reads")]
aligned = re.search('\((.+)%\)', aligned).group(1)
aligned = float(aligned)
mapped = float(values[columns.index("Genome mapped (%)")])
print 'Aligned =', aligned, '%'
print 'Genome mapped =', mapped, '%'
lvl = 0
if aligned < aligned_limit:
    print 'Not enough reads aligned, less than', aligned_limit
    lvl += 1
if mapped < mapped_limit:
    print 'Not enough genome mapped, less than', mapped_limit
    lvl += 2
f.close()
sys.exit(lvl)