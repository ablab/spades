#!/usr/bin/env python3

fin = open("canopy/contigs.in", "r")
fclusters = open("canopy/clusters.out", "r")

fkicked = open("canopy/kicked.log", "w")
fdouble = open("canopy/double.log", "w")

sclusters = set()

while True:
	s = fclusters.readline().strip()
	if s == "":
		break
	cl = s.split()[1]
	if cl in sclusters:
		fdouble.write(cl + "\n")
	else:
		sclusters.add(cl)

while True:
	s = fin.readline()
	if s == "":
		break
	if not s.split()[0] in sclusters:
		fkicked.write(s)


fkicked.close()