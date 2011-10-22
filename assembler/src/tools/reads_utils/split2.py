#!/usr/bin/python -O

import sys

inFile = open(sys.argv[1])
file1 = open(sys.argv[2], 'w')
file2 = open(sys.argv[3], 'w')

id1 = ''
while 1:
	if id1 == '':
		id1 = inFile.readline()

	if not id1:
		break

	l11 = inFile.readline().strip(' \n\t')
	l12 = inFile.readline()

	if l12.startswith('>'):
		id2 = l12
		l12 = '\n'
	else:
		id2 = inFile.readline()

	if (id1.split('/',1)[0] != id2.split('/',1)[0]):
		print 'Non equal'
		break

	file1.write(id1)
	file1.write(l11 + l12)

	l21 = inFile.readline().strip(' \n\t')
	l22 = inFile.readline()

       	if l22.startswith('>'):
               	id1 = l22
               	l22 = '\n'
       	else:
               	id1 = ''

       	file2.write(id2)
       	file2.write(l21	+ l22)		

inFile.close()
file1.close()
file2.close()
