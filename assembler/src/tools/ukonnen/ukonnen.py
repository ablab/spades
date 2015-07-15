#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import itertools

# Changes read format to Ukonnen's scaffolder input format

if len(sys.argv) < 5:
	print 'Prepare reads for Ukonnen\'s scaffolder'
	print 'Usage: python', sys.argv[0], 'INFILE1 INFILE2 OUTFILE1 OUTFILE2'
	exit()

inf1 = open(sys.argv[1]);
inf2 = open(sys.argv[2]);
ouf1 = open(sys.argv[3], 'w');
ouf2 = open(sys.argv[4], 'w');

def reverse_complement(s):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join(map(lambda x: complement[x], s[::-1]))

id = 0
state = 0
for line1, line2 in itertools.izip(inf1, inf2):
	if state == 0: # @...
		assert line1[0] == '@'
		assert line2[0] == '@'
		print >>ouf1, '@' + str(id) + '_R3'
		print >>ouf2, '@' + str(id) + '_F3'
		id += 1
		state = 1
	elif state == 1: # sequence
		assert line1[0] in 'ACGTNacgtn'
		assert line2[0] in 'ACGTNacgtn'
		print >>ouf1, line1.strip()
		print >>ouf2, reverse_complement(line2.strip())
		state = 2
	elif state == 2: # + ...
		assert line1[0] == '+'
		assert line2[0] == '+'
		print >>ouf1, '+'
		print >>ouf2, '+'
		state = 3
	elif state == 3: # quality
		print >>ouf1, line1.strip()
		print >>ouf2, line2.strip()[::-1]
		state = 0
	else:
		assert False
