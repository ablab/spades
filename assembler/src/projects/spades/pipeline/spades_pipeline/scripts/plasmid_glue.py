#!/usr/bin/env python

import os
import sys


def main(args):
    outdir = args[1]
    res = outdir + "/final_contigs.fasta"
    res_f = open(res, "w")
    for file in os.listdir(outdir):
        farr = file.split('.')
        if farr[-1] != "fasta":
            continue
        if farr[-2] != "circular" and farr[-2] != "linears" and farr[-2] != "linearrepeat":
            continue
        
        type = farr[-2]
        arr = farr[-3].split("_")
        if len(arr) < 2:
            continue
        cov = arr[-1]
        if len(cov) > 4: 
            continue

  #  for line in open(os.path.join(dir,file), "r"):
        print (file)
        for line in open(os.path.join(outdir,file), "r"):
            line = line.strip()
            if len(line) > 0 and line[0] == ">":
                line += "_cutoff_" + cov+ "_type_" + type
            res_f.write(line+ "\n")
    res_f.close()
    if os.path.getsize(res) != 0:        
        scaff = outdir + "/scaffolds.fasta"   
        from shutil import copyfile
        copyfile(res, scaff)
    else:
        os.remove(res)
if __name__ == "__main__":
    main(sys.argv)

