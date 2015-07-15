#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#
# This script removes whitespaces from headers in FASTA-file 
#

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
