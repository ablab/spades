#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
file1 = open(sys.argv[2], 'w')
file2 = open(sys.argv[3], 'w')

format = 'fasta'
lines = 1
if len(sys.argv) >= 5:
	format = sys.argv[4]

if format == 'fastq':
	lines = 2
	
done = 2
while 1:
	if done == 2:
		line = inFile.readline()

	if not line:
		break
	
	done = 0
	l1 = ""
	for i in range(0,lines):
		l1 = inFile.readline()
		if not l1:
			line = l1
			break

		if l1.startswith('@') or l1.startswith('>'):
			line = l1
			break
		done = 1


	if done == 1:
		for i in range(0,lines):
			l = inFile.readline()
			if not l:
				line = l
				break

			if l.startswith('@') or l.startswith('>'):
				line = l
				break

			file1.write(line.strip(' \t\n') + '/1\n') 
			file1.write(l1)

			file2.write(line.strip(' \t\n') + '/2\n') 
			file2.write(l)
			done = 2


		
inFile.close()
file1.close()
file2.close()
