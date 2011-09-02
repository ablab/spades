#!/usr/bin/python

import sys
import re

if len(sys.argv) < 2:
  print 'Usage:', sys.argv[0], 'FASTA_FILE'
  exit(1)

fastafile = sys.argv[1]
for line in open(fastafile):
  if line[0] == '>':
    print re.sub(r'\s', '', line)
  else:
    print line,
