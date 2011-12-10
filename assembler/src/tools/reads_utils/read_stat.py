#!/usr/bin/python -O

#Calculate statistics for reads by reference file

import argparse
import sys
import os

argParser = argparse.ArgumentParser(description='Calculate reads statistic using bowtie aligner')
argParser.add_argument('bowtie', metavar='bowtie_path', type=str, help='path to bowtie aligner')
argParser.add_argument('index', metavar='genome_index', type=str, help='genome bowtie index')
argParser.add_argument('length', metavar='genome_length', type=int, help='genome length')
argParser.add_argument('--chrs', type=str, help='chromosomes map file')
argParser.add_argument('--mode', type=str, default="fr",  choices=["fr", "rf", "ff", "s"], help='alignment mode')
argParser.add_argument('-p', type=int, default=8, help='threads used by bowtie, default 8')
argParser.add_argument('--fst', type=str, help='1st reads file')
argParser.add_argument('--snd', type=str, help='2nd reads file, can be ommited if only one single read file is used')
argParser.add_argument('--minis', type=int, default=0, help='minimal insert size')
argParser.add_argument('--maxis', type=int, default=500, help='maximal insert size')
argParser.add_argument('--separate_alignment', action='store_true', help='separately align paired reads and than calculate paired statistics, can be more useful')
argParser.add_argument('--use_multiple', action='store_true', help='use multiple aligned reads, does not work if --separate_alignment is used')
argParser.add_argument('--leave_tmp', action='store_true', help='do not clear temperary file')

args = argParser.parse_args()

print(args.use_multiple)

suppress_muliple = " -m 1 "
if args.use_multiple:
	suppress_muliple = " "

fName, ext = os.path.splitext(args.fst)

format = 'q'
if ext == '.fasta':
	format = 'f'

outp = ""

if args.mode != "s" and not args.separate_alignment:
	outp = args.fst
	os.system(args.bowtie + ' -c -' + format + ' -p ' + str(args.p) + suppress_muliple + ' --' + args.mode + ' -I ' + args.minins + ' -X ' + args.maxins + ' --suppress 6,7,8 ' 
			+ args.index + ' -1 ' + args.fst + ' -2 ' + args.snd + ' > ' + outp + '.log 2> ' + outp + '.stat') 

elif args.separate_alignment:
	os.system(args.bowtie + ' -c -' + format + ' -p ' + str(args.p) + ' --best ' + ' --suppress 6,7,8 ' 
			+ args.index + ' ' + args.fst + ' > ' + args.fst + '.log 2> ' + args.fst + '.stat')
	os.system(args.bowtie + ' -c -' + format + ' -p ' + str(args.p) + ' --best ' + ' --suppress 6,7,8 ' 
			+ args.index + ' ' +args.snd + ' > ' + args.snd + '.log 2> ' + args.snd + '.stat')

else:	
	outp = args.fst + '.tmp' + ext
	os.system('cp ' + args.fst + ' > ' + outp)
	os.system('cat ' + args.snd + ' >> ' + outp)
	os.system(args.bowtie + ' -c -' + format + ' -p ' + str(args.p) + suppress_muliple + ' --suppress 6,7,8 ' 
			+ args.index + ' ' + outp + ' > ' + outp + '.log 2> ' + outp + '.stat')


#bowtie -c -q -p 8 -m 1 --rf -I 0 -X 2000 --suppress 6,7,8 index/SAureus -1 JCVI_control_MDA_Saureus_10.1.fastq  -2 JCVI_control_MDA_Saureus_10.2.fastq > sa_rf.log 2> sa_rf.err &
	

