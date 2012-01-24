#!/usr/bin/python

import sys
import os
import shutil
import re
import subprocess

f = open(sys.argv[1], 'r')
columns = map(lambda s: s.strip(), f.readline().split('\t'))
values = map(lambda s: s.strip(), f.readline().split('\t'))
n50 = int(values[columns.index("N50")])
mis = int(values[columns.index("Misassemblies")])
print 'n50 = ', n50
print 'missasembled contigs = ', mis
if n50 < 75000:
    print('N50 is too small')
    sys.exit(1)
if (mis > 2):
    print('too many missassemblies')
    sys.exit(1)
f.close()