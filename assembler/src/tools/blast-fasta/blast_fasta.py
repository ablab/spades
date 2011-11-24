#!/usr/bin/python

import sys

# BLASTN all items in FASTA file

if len(sys.argv) < 2:
	print 'Usage: python', sys.argv[0], ' FILE_NUCL.FASTA'
	exit()

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import re

records = list(SeqIO.parse(sys.argv[1], 'fasta', generic_dna)) # all FASTA to memory (works on bacteria well)

l = len(records)
for i, record in enumerate(records):
	print str(i+1) + '/' + str(l) + ':', record.id
	result_handle = NCBIWWW.qblast("blastn", "nr", record.seq)
	blast_record = NCBIXML.read(result_handle)
	E_VALUE_THRESH = 0.04
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < E_VALUE_THRESH:
				print '****Alignment****'
				print 'sequence:', alignment.title
				print 'length:', alignment.length
				print 'e value:', hsp.expect
				print hsp.query
				print hsp.match
				print hsp.sbjct
	print

