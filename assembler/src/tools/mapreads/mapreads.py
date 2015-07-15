#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import sys
import os

BOWTIE_DIR = '../../../ext/tools/bowtie-0.12.7/'

if len(sys.argv) < 6:
	print 'Usage:', sys.argv[0], ' FASTA_REFERENCE  FASTQ_READS_1  FASTQ_READS_2  MIN_INSERT_SIZE  MAX_INSERT_SIZE'
	print 'Maps short paired reads to reference (e.g. contigs), produces sam file'
	exit(1)

reference = sys.argv[1]
reads1, reads2 = sys.argv[2:4]
min_insert_size = sys.argv[4]
max_insert_size = sys.argv[5]

os.system(BOWTIE_DIR + '/bowtie-build ' + reference + ' ./ebwt_outfile_base')
os.system(BOWTIE_DIR + '/bowtie ./ebwt_outfile_base -1 ' + reads1 + ' -2 ' + reads2 + 
           ' --sam ' + ' --minins ' + min_insert_size + ' --maxins ' + max_insert_size)
os.system('rm ./ebwt_outfile_base*')
