#!/usr/bin/env python3

import os
import sys

CONTIGS_MPL_DIR = "contigs-mpl"
CANOPY = "/home/toxa31/libs/canopy/cc.bin"
OUTPUT_DIR = "canopy"
INPUT_FILE = "contigs.in"

files = [f[:-3] for f in os.listdir(CONTIGS_MPL_DIR) if os.path.isfile(os.path.join(CONTIGS_MPL_DIR, f)) and f[-3:] == ".id"]

if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)
else:
	os.system("rm -rf " + OUTPUT_DIR + "/*")

inp = open(os.path.join(OUTPUT_DIR, INPUT_FILE), "w")

for f in files:
	ctg_id = open(os.path.join(CONTIGS_MPL_DIR, f + ".id"), "r")
	ctg_mpl = open(os.path.join(CONTIGS_MPL_DIR, f + ".mpl"), "r")
	while True:
		cid = ctg_id.readline().strip()
		cmpl = ctg_mpl.readline().strip()
		if ((cid == "") != (cmpl == "")):
			print("Unexpected end of file")
			sys.exit(1)
		if (cid == ""):
			break

		inp.write(f + "-" + cid + " " + cmpl + "\n")

	ctg_id.close()
	ctg_mpl.close()

inp.close()
