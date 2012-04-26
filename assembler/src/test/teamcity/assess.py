#!/usr/bin/python

############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
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
mis = int(values[columns.index("Misassemblies")])
print 'N50 =', n50
print 'Misasembled contigs =', mis
if n50 < n50_limit:
    print 'N50 is too small: less than', n50_limit
    sys.exit(1)
if mis > mis_limit:
    print 'Too many misassemblies: more than', mis_limit
    sys.exit(1)
f.close()
