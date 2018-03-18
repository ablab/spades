#!/usr/bin/env python

import os
import sys
cur_d = os.path.dirname(os.path.realpath(__file__))

bin = os.path.join(cur_d, "bin/spades")

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " config file "
dir = ""
K = ""
os.system (bin  + " " + sys.argv[1] + " " + os.path.join(cur_d, "configs/debruijn/mda_mode.info") + " " +  os.path.join(cur_d, "configs/debruijn/meta_mode.info") + " " +  os.path.join(cur_d, "configs/debruijn/plasmid_mode.info"))
for line in open(sys.argv[1], 'r'):
    arr = line.split()
    if len(arr) >= 2:
        if arr[0] == "K":
            K = "K" + arr[1]
        if arr[0] == "output_base":
            dir = arr[1]
dir = dir + "/" + K
res = dir + "/final_contigs_cat.fasta"
res_f = open(res, "w")
for file in os.listdir(dir):
    farr = file.split('.')
    if farr[-1] != "fasta":
        continue
    if farr[-2] != "circular":
        continue
    cov = farr[-3].split("_")[-1]
    for line in open(os.path.join(dir,file), "r"):
        line = line.strip()
        if len(line) > 0 and line[0] == ">":
            line += "_cutoff_" + cov
        res_f.write(line+ "\n")
