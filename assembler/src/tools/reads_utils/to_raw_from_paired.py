#!/usr/bin/python -O

import sys
import os

inFileName = sys.argv[1]
inFile= open(inFileName, 'r')

fName, ext = os.path.splitext(inFileName)
frFile = open(fName + "_FR" + ext, "w")
rfFile = open(fName + "_RF" + ext, "w")
ffFile = open(fName + "_FF" + ext, "w")
 
delim = '/'

ffC = 0
frC = 0
rfC = 0

prevLine = ""
for line in inFile:
	if (prevLine == ""):
		prevLine = line
		continue

	if (line.split(delim, 1)[0] == prevLine.split(delim, 1)[0]):
		l1 = line.split('\t', 5)
		l2 = prevLine.split('\t', 5)

		or2 = l1[1]
		or1 = l2[1]
		pos2 = int(l1[3])
		pos1 = int(l2[3])
               	len2 = len(l1[4])
               	len1 = len(l2[4])

		if (or1 == or2):
			ffC += 1
			if (pos1 < pos2):
				ffFile.write(str(pos1) + ' ' + str(len1) + '\n' + str(pos2) + ' ' + str(len2) + '\n')
			else:
				ffFile.write(str(pos2) + ' ' + str(len2) + '\n' + str(pos1) + ' ' + str(len1) + '\n')

		elif (pos1 <= pos2 and or1 == '+'):
			frC += 1
			frFile.write(str(pos1) + ' ' + str(len1) + '\n' + str(pos2) + ' ' + str(len2) + '\n')

		elif (pos2 < pos1 and or2 == '+'):
                        frC += 1
                        frFile.write(str(pos2) + ' ' + str(len2) + '\n' + str(pos1) + ' ' + str(len1) + '\n')

                elif (pos1 <= pos2 and or1 == '-'):
                        rfC += 1
                        rfFile.write(str(pos1) + ' ' + str(len1) + '\n' + str(pos2) + ' ' + str(len2) + '\n')

                elif (pos2 < pos1 and or2 == '-'):
                        rfC += 1
                        rfFile.write(str(pos2) + ' ' + str(len2) + '\n' + str(pos1) + ' ' + str(len1) + '\n')

		else:
			print('Something wrong: '+str(pos1)+' '+str(len1)+' '+or1+' '+str(pos2)+' '+str(len2)+' '+or2+'\n')

	else:
		print("Non-equal pairs\n")
		print(prevLine.split(delim, 1)[0])
		print(line.split(delim, 1)[0])

	prevLine = ""


print('FR: '+str(frC)+'\nRF: '+str(rfC)+'\nFF: '+str(ffC)+'\n')

inFile.close()
frFile.close()
rfFile.close()
ffFile.close()
