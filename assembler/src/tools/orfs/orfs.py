#!/usr/bin/python
import sys
import itertools

sys.path.append('../quality/libs')
import fastaparser

if len(sys.argv) < 3:
	print 'Counts bacterial ORFs (trans table 11)'
	print 'Usage: python', sys.argv[0], ' FASTA_FILE ORF_MIN [ORF_MAX]'
	print 'ORFs length are in codons (3bp), including start and stop codons.'
	exit()

# http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
# 11. The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
Starts = '---M---------------X------------XXXM---------------M------------' # Ms are starts, Xs are deprecated (rare) starts
Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

fasta_filename = sys.argv[1]
orf_min = int(sys.argv[2])
orf_max = int(sys.argv[3]) if len(sys.argv) >= 4 else 1e3000

# precalc of start/stop codons based on transl_table
starts = set()
stops = set()
for i in xrange(64):
	if AAs[i] == '*':
		stops.add((Base1[i], Base2[i], Base3[i]))
	if Starts[i] == 'M':
		starts.add((Base1[i], Base2[i], Base3[i]))

def find_ORFs(genome, start, min, max):
	gene = False # no gene
	l = 0 # length of ORF
	ls = []
	for a, b, c in itertools.izip(genome[start::3], genome[start+1::3], genome[start+2::3]):
		l += 1
		if not gene and (a, b, c) in starts:
			gene = True # gene
			l = 1
		elif gene and (a, b, c) in stops:
			gene = False
			if min <= l <= max:
				ls.append(l)
	return ls

def reverse_complement(s):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join(map(lambda x: complement[x], s[::-1]))

fasta = fastaparser.read_fasta(fasta_filename)
cnt = 0
for name, seq in fasta:
	seq = seq.upper()
	cnt += len(find_ORFs(seq, 0, orf_min, orf_max))
	cnt += len(find_ORFs(seq, 1, orf_min, orf_max))
	cnt += len(find_ORFs(seq, 2, orf_min, orf_max))
	rc_seq = reverse_complement(seq)
	cnt += len(find_ORFs(rc_seq, 0, orf_min, orf_max))
	cnt += len(find_ORFs(rc_seq, 1, orf_min, orf_max))
	cnt += len(find_ORFs(rc_seq, 2, orf_min, orf_max))

print 'ORFs between', orf_min, 'and', orf_max, 'in', fasta_filename,' = ', cnt