#! /usr/bin/env python
import sys
import argparse
#import os
import subprocess
from joblib import Parallel, delayed
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
#from Bio import SeqIO
from glob import glob


def parse_args(args):
    ###### Command Line Argument Parser
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument('--b', help='Path to blast database')
    parser.add_argument('--f', help='Path to folder with SPAdes outputs')   
    return parser.parse_args()


def main():
	args = parse_args(sys.argv[1:])

# get the files

	scaffold_files = glob(str(args.f)+"/*/scaffolds.fasta")
	print (scaffold_files)

# for each scaffolds file - blast, save as xml
	blastn_args_list = []
	for i in scaffold_files:
		blastn_args_list.append(["blastn","-query", i.strip(), "-db", args.b, "-evalue", "0.001", "-outfmt", "5", "-out",str(i.split("/")[-2])+".xml", "-num_alignments","20"])
		print (blastn_args_list[-1])
	Parallel(n_jobs=30)(delayed(subprocess.call) (args) for args in blastn_args_list)
		
'''
		blastn_cline = NcbiblastnCommandline(query=i, db=args.b, evalue=0.001, outfmt=5, out=str(i.split("/")[-2])+".xml", num_alignments=20)
		print (blastn_cline)
		blastn_cline
		stdout, stderr = blastn_cline()
'''




if __name__ == '__main__':
    main()


