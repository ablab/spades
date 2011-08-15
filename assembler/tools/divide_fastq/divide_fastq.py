import sys

if len(sys.argv) < 2:
	print 'Usage:', sys.argv[0], ' INFILE'
	exit(0)
	
infilename = sys.argv[1]
outfilename1 = infilename + '_1'
outfile1 = open(outfilename1, 'w')
outfilename2 = infilename + '_2'
outfile2 = open(outfilename2, 'w')
cnt = 0
for line in open(filename):
	if cnt % 8 < 4:
		print >>outfile1, line
	else:
		print >>outfile2, line
	cnt += 1
