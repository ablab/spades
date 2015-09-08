#!/usr/bin/env python3

import os

COUNT_SAMPLES = 4
REFS_DIR = "refs"
OUT_DIR = "gen-samples"
WGSIM_PATH = "/home/toxa31/libs/wgsim/wgsim"

files = [f for f in os.listdir(REFS_DIR) if os.path.isfile(os.path.join(REFS_DIR, f))][:COUNT_SAMPLES]

if not os.path.exists(OUT_DIR):
	os.makedirs(OUT_DIR)
else:
	os.system("rm -rf " + OUT_DIR + "/*")	


for (i, f) in enumerate(files):
	print("Generating reads for {} [{}/{}]".format(f, i + 1, len(files)))	
	fname = os.path.splitext(f)[0]
	os.system("{wgsim} -N 100000 -r 0.01 -1 100 -2 100 -S31 -d0 -e0 {refs}/{ref} {out}/{read}.r1.fastq {out}/{read}.r2.fastq > /dev/null 2>&1"
		.format(wgsim=WGSIM_PATH, refs=REFS_DIR, out=OUT_DIR, ref=f, read=fname))