#!/usr/bin/env python
from __future__ import print_function

import argparse
import os.path

argparser = argparse.ArgumentParser(description="Make csv summary from clusters")
argparser.add_argument("--output", "-o", type=str, help="Output csv")
argparser.add_argument("input", type=str, help="Directory with binning info")

from os import listdir
from os.path import isfile, join


if __name__ == "__main__":
	args = argparser.parse_args()
	mypath = "/".join(args.input.split("/")[:-1])
	fastafiles = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f)) and f.endswith("fasta")]
	res = ""
	for f in fastafiles:
	    fin = open(f, "r")
	    cluster_num = int(f.split(".")[-2])
	    for l in fin.readlines():
	        if l.startswith(">"):
	            contig = l.strip()[1:]
	            res += contig + "," + str(cluster_num) + "\n"
	    fin.close()
	fout = open(args.output, "w")
	fout.write(res)
	fout.close()
