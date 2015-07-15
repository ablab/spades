#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


#Calculate statistics for reads by reference file

import argparse
import sys
import os

argParser = argparse.ArgumentParser(description='Calculate reads statistic using bowtie aligner')
argParser.add_argument('bowtie', metavar='bowtie_path', type=str, help='path to bowtie aligner')
argParser.add_argument('genome', metavar='genome', type=str, help='genome reference or contigs')
argParser.add_argument('--genome_length', metavar='genome_length', type=int, default=0, help='genome length (for counting coverage)')
#argParser.add_argument('--chrs', type=str, help='chromosomes map file')
argParser.add_argument('--mode', type=str, default="fr",  choices=["fr", "rf", "ff", "s"], help='alignment mode')
argParser.add_argument('-p', type=int, default=8, help='threads used by bowtie, default 8')
argParser.add_argument('--fst', type=str, help='1st reads file')
argParser.add_argument('--snd', type=str, help='2nd reads file, can be ommited if only one single read file is used')
argParser.add_argument('--minis', type=int, default=0, help='minimal insert size')
argParser.add_argument('--maxis', type=int, default=500, help='maximal insert size')
#argParser.add_argument('--separate_alignment', action='store_true', help='separately align paired reads and than calculate paired statistics, can be more useful')
#argParser.add_argument('--use_multiple', action='store_true', help='use multiple aligned reads, does not work if --separate_alignment is used')
argParser.add_argument('--leave_tmp', action='store_true', help='do not clear temperary file')

args = argParser.parse_args()

outdir = 'tmp__'
os.system('mkdir -p ' + outdir + '/index')
os.system(args.bowtie + 'bowtie-build ' + args.genome + ' ' + outdir + '/index/ref' + ' > /dev/null')

fName, ext = os.path.splitext(args.fst)

format = 'q'
if ext == '.fasta':
	format = 'f'

outp =  outdir + '/' + args.fst + '.tmp' + ext

if args.mode != "s":
	os.system(args.bowtie + 'bowtie -c -' + format + ' -p ' + str(args.p) + ' --best --' + args.mode + ' -I ' + args.minis + ' -X ' + args.maxis + ' --suppress 6,7,8 ' 
			+ outdir + '/index/ref' + ' -1 ' + args.fst + ' -2 ' + args.snd + ' > ' + outp + '.log 2> ' + outp + '.stat') 

	try:
		#Redirecting output
		#bufferString = StringIO()
		#sys.stdout = bufferString
		#sys.stderr = bufferString
		#Executing
		execfile("stat/raw.py " + outp + '.log ' + outp + '.raw' )
		execfile("stat/is.py " + outp + '.raw ' + outp + '.is ' + str(args.minis) + ' ' + str(args.maxis) )
		execfile("stat/mean.py " + outp + '.is ')
		execfile("stat/coverage.py " + outp + '.raw ' + outp + '.cov ' + str(genome_length) + ' 1000 ')

	

	except IOError:
		pass

	finally:
		#Returning output to its standart values
		sys.stdout = sys.__stdout__
		sys.stderr = sys.__stderr__


else:	
	os.system('cp ' + args.fst + ' ' + outp)
	if args.snd is not None:
		os.system('cat ' + args.snd + ' >> ' + outp)

	os.system(args.bowtie + 'bowtie -c -' + format + ' -p ' + str(args.p) + suppress_muliple + ' --suppress 6,7,8 ' 
			+ args.index + ' ' + outp + ' > ' + outp + '.log 2> ' + outp + '.stat')

	try:
		execfile("stat/raw.py " + outp + '.log ' + outp + '.raw' )
		execfile("stat/coverage.py " + outp + '.raw ' + outp + '.cov ' + str(genome_length) + ' 1000 ')

	except IOError:
		pass

	finally:
		#Returning output to its standart values
		sys.stdout = sys.__stdout__
		sys.stderr = sys.__stderr__


if not args.leave_tmp:
	os.system('rm -r tmp__')
	

