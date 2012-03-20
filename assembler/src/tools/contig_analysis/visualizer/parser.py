#!/usr/bin/python -O

import sys
import os

if len(sys.argv) != 2:
	print("Usage: " + sys.argv[0] + " <plantagota output>")	
	sys.exit()

inFileName = sys.argv[1]
inFile = open(inFileName)

outFileName, ext = os.path.splitext(inFileName)
outFileName += ".txt"

print (outFileName)

outFile = open(outFileName, "w")

for line in inFile:
	parse = line.strip().split(' ')
	if parse[0] == "Align":
		forward = " +"
#		print(parse[2] + " " + parse[3] + " " + parse[4])
		if (int(parse[2]) - int(parse[3]))*(int(parse[5]) - int(parse[6])) < 0:
			forward = " -"
		outFile.write(parse[2] + " " + parse[3] + " " + parse[4] + forward + " " + str(min(int(parse[5]),int(parse[6]))) + "\n")

outFile.close()
inFile.close()
		
	
