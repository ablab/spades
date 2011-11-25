#!/usr/bin/python -O

#Calculate statistics for reads by reference file

import argparse
import sys
import os

argParser = argparse.ArgumentParser(description='Calculate read statistic using bowtie aligner')
argParser.add_argument('bowtie', metavar='bowtie_path', type=str, help='path to bowtie aligner')
argParser.add_argument('index', metavar='genome_index', type=str, help='genome bowtie index')
argParser.add_argument('length', metavar='genome_length', type=int, help='genome length')
argParser.add_argument('--mode', type=str, default="fr",  choices=["fr", "rf", "ff", "s"], help='alignment mode')
argParser.add_argument('--fst', type=str, help='1st reads file')
argParser.add_argument('--snd', type=str, help='2nd reads file, can be ommited if only one single read file is used')
argParser.add_argument('--minis', type=int, default=0, help='minimal inser size')
argParser.add_argument('--maxis', type=int, default=500, help='maximal insert size')
argParser.add_argument('--separate_alignment', action='store_true', help='separately align paired reads and than calculat paired statistics, can be more useful')
argParser.add_argument('--use_multiple', action='store_true', help='use multiple aligned reads, does not work if --separate_alignment is used')
argParser.add_argument('--leave_tmp', action='store_true', help='do not clear temperary file')

args = argParser.parse_args()

print(args.use_multiple)

if args.mode != s and not args.separate_alignment:
	os.system(args.bowtie + ' ' + 



	

