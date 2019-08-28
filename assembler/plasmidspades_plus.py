#!/usr/bin/env python

import os, errno
import sys
import argparse
import collections
from math import log
from math import exp
import csv
import operator



def parse_args(args):
###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Wrapper script for plasmidSPAdes accomplished with plasmidVerify verification")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument('-d', required = True, help='Input dataset description file(in yaml format, see SPAdes manual for details)')
    parser.add_argument('-o', required = True, help='Output directory')
    parser.add_argument('--pv', required = True, help='Path to plasmidVerify script')
    parser.add_argument('--db', help='Path to BLAST db')
    parser.add_argument('--hmm', required = True, help='Path to Pfam-A HMM database')    
    return parser.parse_args()



args = parse_args(sys.argv[1:])
current_path = os.path.dirname(os.path.realpath(__file__))
spades = os.path.join(current_path, "spades.py")
plasmidverify = os.path.join(args.pv, "plasmidverify.py")

spades_string = spades + " --plasmid --dataset " + args.d + " -o " + args.o 
os.system(spades_string)
plasmidverify_output = os.path.join (args.o, "plasmidverify")
scaffolds = os.path.join(args.o, "scaffolds.fasta")
plasmidverify_string = "python " + plasmidverify + " -o " + plasmidverify_output + " --hmm " + args.hmm + " -f " + scaffolds
if args.db:
    plasmidverify_string += " --db " + args.db
os.system(plasmidverify_string)

