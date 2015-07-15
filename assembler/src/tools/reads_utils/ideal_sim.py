#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# makes single reads of length RL
RL = 220

f = open('MG1655-K12.fasta')
f.readline()
s = ''
for line in f:
  s += line.strip()
l = len(s)
s += s[:RL] # make it circular
for i in xrange(l):
  print '@read%d' % i 
  print s[i:i+RL] 
  print '+'
  print 'B' * RL
