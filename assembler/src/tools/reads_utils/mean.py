#!/usr/bin/python -O

import sys
import math

inFile = open(sys.argv[1])

ttl = 0.0
cnt = 0.0

while (1):
	line = inFile.readline()

	if not line:
		break

	ttl += int(line.split(' ')[0]) * int(line.split(' ')[1])
	cnt += int(line.split(' ')[1])

mean = ttl/cnt

inFile.close()
inFile = open(sys.argv[1])

dev = 0.0
while (1):
       	line = inFile.readline()

       	if not line:
               	break

       	dev += (int(line.split(' ')[0]) - mean) * (int(line.split(' ')[0]) - mean) * int(line.split(' ')[1])

print("Mean: " + str(mean) + "\nDeviation: " + str(math.sqrt(dev / cnt)))

inFile.close()
