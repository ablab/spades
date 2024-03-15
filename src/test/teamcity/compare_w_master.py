#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################



import sys
import os
import glob



def find_in_dir(contig_dir, prefix):
    master_ctgs = sorted(glob.glob(os.path.join(contig_dir, "*_" + prefix + "_master.fasta")))
    other_ctgs = sorted(glob.glob(os.path.join(contig_dir, "*_" + prefix + ".fasta")))

    other_ctg = other_ctgs[-1]
    if other_ctg.find("latest") != -1:
        other_ctg = other_ctgs[-2]
    if len(master_ctgs) > 0 and len(other_ctgs) > 0:
        return [master_ctgs[-1], other_ctg]
    return []

# Run QUAST for a set of contigs
def run_quast(contigs, quast_output_dir, ref, opts):
    if not reduce(lambda x, y: os.path.exists(y) and x, contigs, True):
        log.warn("No contigs were found in " + output_dir)
        return 8

    cmd = "/home/teamcity/quast_latest/quast.py -R " + ref + " " + opts + " "
    print('Running ' + cmd + ' on ' + ','.join(contigs))
    quast_cmd = cmd + " -o " + quast_output_dir + " " + " ".join(contigs) + " > /dev/null"
    
    ecode = os.system(quast_cmd)
    if ecode != 0:
        print(quast_cmd + " finished abnormally with exit code " + str(ecode))
        return 9

    return 0


def run_all_quasts(contig_dir, ref, quast_output_dir):
    i = 1
    for prefix, opts in [("contigs", ""), ("scaffolds", ""), ("scaffolds", " -s ")]:
        contigs = find_in_dir(contig_dir, prefix)
        if len(contigs) == 2:
            run_quast(contigs, os.path.join(quast_output_dir, str(i) + "_" + prefix), ref, opts)
        else:
            print("Skipping contigs " + " ".join(contigs))
        i += 1


if len(sys.argv) != 4:
    print("Run: " + sys.argv[0] + " <contig dir> <refrence> <output>")
    exit(1)

run_all_quasts(sys.argv[1], sys.argv[2], sys.argv[3])

            
