#!/usr/bin/env python3

import random

ctg = open("canopy/contigs.in", "r")
ctr = open("canopy/clusters.out", "r")

out = open("canopy/points_matrix.csv", "w")

ctg_to_ctr = dict()

while True:
	s = ctr.readline().strip()
	if (s == ""):
		break
	a = s.split()
	ctr_id = a[0][3:]

	if (random.randint(1, 25) == 1):
		ctg_to_ctr[a[1]] = ctr_id

while True:
	s = ctg.readline().strip()
	if s == "":
		break

	a = s.split()
	if (a[0] in ctg_to_ctr):
		out.write(ctg_to_ctr[a[0]])
		for x in a[1:]:
			out.write("," + x)

		out.write("\n")

out.close()