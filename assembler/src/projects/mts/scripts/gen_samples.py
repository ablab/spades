#!/usr/bin/env python3

import os
import random

COUNT_SAMPLES = 5
COUNT_REFS = 2
REFS_DIR = "/Sid/eplekhanova/refs"
OUT_DIR = "/Sid/eplekhanova/samples"
WGSIM_PATH = "wgsim"

files = [f for f in sorted(os.listdir(REFS_DIR)) if os.path.isfile(os.path.join(REFS_DIR, f))][:COUNT_REFS]

if os.path.exists(OUT_DIR):
	os.system("rm -rf " + OUT_DIR + "/*")


for i in range(1, COUNT_SAMPLES + 1):
	print("Generating reads for sample {}/{}".format(i + 1, COUNT_SAMPLES))

	os.system("mkdir -p {}/sample{}".format(OUT_DIR, i))

	for f in files:
		fname = os.path.splitext(f)[0]
		os.system("{wgsim} -N {reads} -r 0.01 -1 100 -2 100 -S{seed} -d 300 -e0 {refs}/{ref} {out}/{read}.tmp.r1.fastq {out}/{read}.tmp.r2.fastq > /dev/null 2>&1"
			.format(wgsim=WGSIM_PATH, refs=REFS_DIR, out=OUT_DIR, ref=f, read=fname, seed=i, reads=int((10 ** 6) * random.uniform(0.25, 2.0))))

	os.system("cat {out}/*.tmp.r1.fastq >> {out}/sample{ind}/r1.fastq".format(out=OUT_DIR, ind=i))
	os.system("cat {out}/*.tmp.r2.fastq >> {out}/sample{ind}/r2.fastq".format(out=OUT_DIR, ind=i))
	os.system("rm -f {out}/*.tmp.*.fastq".format(out=OUT_DIR))
