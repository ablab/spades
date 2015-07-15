############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

def get_lengths_from_fastafile(filename):
	lengths = []
	l = 0
	for line in open(filename):
		if (line[0] == '>'):
			if l != 0: # not first sequence in fasta
				lengths.append(l)
				l = 0
		else:
			l += len(line.strip())
	lengths.append(l)
	return lengths
