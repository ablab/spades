#!/usr/bin/env python3

import os
import subprocess

SAMPLES_DIR = "gen-samples"
SPADES = "/home/toxa31/work/algorithmic-biology/assembler/spades.py"
OUTPUT_DIR = "contigs-mpl"

files = [f[:-9] for f in os.listdir(SAMPLES_DIR) if os.path.isfile(os.path.join(SAMPLES_DIR, f))]
seen = set()
for f in files:
	seen.add(f)
files = list(seen)

if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)
else:
	os.system("rm -rf " + OUTPUT_DIR + "/*")

for f in files:
	cmd = "{spades} --only-assembler -1 {dir}/{read}.r1.fastq -2 {dir}/{read}.r2.fastq -o /home/toxa31/work/assembled".format(
		spades=SPADES, read=f, dir=SAMPLES_DIR)
	if subprocess.call(cmd, shell=True) != 0:
		print("")
		print("Finished with error!")
		break


	os.system("mv contigs.id " + OUTPUT_DIR + "/" + f + ".ctg.id")
	os.system("mv contigs.mpl " + OUTPUT_DIR + "/" + f + ".ctg.mpl")
