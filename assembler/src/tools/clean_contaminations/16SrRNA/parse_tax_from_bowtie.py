#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


# Before run: 
# bowtie -p 6 --all --best --chunkmbs 1024 SSURef_108_tax_silva_trunc_DNA ecoli_mda_lane1.fastq > ecoli_mda_lane1_bowtie_allbest.sam
# bowtie -p 6 -v 1 --all --best --chunkmbs 1024 SSURef_108_tax_silva_trunc_DNA ecoli_mda_lane1.fastq > ecoli_mda_lane1_bowtie_allbest.sam
# bowtie -p 6 -k 50 --best --chunkmbs 1024 SSURef_108_tax_silva_trunc_DNA ecoli_mda_lane1.fastq > ecoli_mda_lane1_bowtie_k50best.sam

# Gets taxonomy dominants from bowtie/bwa sam-file

import sys

MIN_FREQUENCY = 0.90 # from 0 to 1

# matches = []
taxes = {}
total = 0
for samline in sys.stdin:
		line = samline.split('\t')
		assert line[1] == '+' or line[1] == '-', line # correct sam file
		line = line[2] # fasta id for match
		line = line[line.find('_')+1:] # throw out read id (only taxonomy left)
		# matches.append(line)
		line = line.split(';')
		total += 1
		key = ''
		for tax in line:
			key += tax
			if key in taxes:
				taxes[key] += 1
			else:
				taxes[key] = 1
			key += ';'

taxes_ls = zip(taxes.values(), taxes.keys())
taxes_ls.sort(key = lambda x: (x[0], x[1] + chr(255)), reverse = True)

#top = 10
#for tax in taxes_ls[:top]:
#	print tax[0], tax[1]

for tax in taxes_ls:
    print tax

print

choice = ''
depth = 1
for tax in taxes_ls:
	name = tax[1]
	ls = name.split(';')
	if len(ls) == depth:
		frequency = tax[0] * 1.0 / total
		print '%6.2f%%:' % (frequency * 100), name
		if frequency >= MIN_FREQUENCY:
			choice = name
		depth += 1
print
print 'Choice (with >= %.2f%%):' % (MIN_FREQUENCY * 100), (choice if choice else 'NO')

