############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import os

if len(sys.argv) < 2:
	print 'Usage:', sys.argv[0], ' INFILE'
	exit(0)
	
infilename = sys.argv[1]
outfilename1 = os.path.basename(infilename) + '_1'
outfile1 = open(outfilename1, 'w')
outfilename2 = os.path.basename(infilename) + '_2'
outfile2 = open(outfilename2, 'w')
cnt = 0
for line in open(infilename):
	if cnt % 8 < 4:
		print >>outfile1, line,
	else:
		print >>outfile2, line,
	cnt += 1
