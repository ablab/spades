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

AAs    = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
Starts = '---M---------------M------------MMMM---------------M------------'
Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'

fasta_filename = sys.argv[1]
orf_min = int(sys.argv[2])
orf_max = int(sys.argv[3]) if len(sys.argv) >= 4 else 1e3000
fasta = fastaparser.read_fasta(fasta_filename)

def find_ORFs(genome):
	gene = False # no gene
	l = 0 # length of ORF
	ls = []
	for a, b, c in itertools.izip(genome[::3], genome[1::3], genome[2::3]):
		l += 1
		for i in xrange(64):
			if (Base1[i] == a) and (Base2[i] == b) and (Base3[i] == c):
				if Starts[i] == 'M' and not gene:
					gene = True # gene
					l = 1
				elif AAs[i] == '*' and gene:
					gene = False
					ls.append(l)
				break
	return ls

def cnt_ORFs(genome, min, max):
	orfs = find_ORFs(genome)
	return len(filter(lambda x: min <= x <= max, orfs))

def reverse_complement(s):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
	return ''.join(map(lambda x: complement[x], s[::-1]))

cnt = 0
for name, seq in fasta:
	cnt += cnt_ORFs(seq, orf_min, orf_max)
	cnt += cnt_ORFs(seq[1:], orf_min, orf_max)
	cnt += cnt_ORFs(seq[2:], orf_min, orf_max)
	rc_seq = reverse_complement(seq)
	cnt += cnt_ORFs(rc_seq, orf_min, orf_max)
	cnt += cnt_ORFs(rc_seq[1:], orf_min, orf_max)
	cnt += cnt_ORFs(rc_seq[2:], orf_min, orf_max)


print 'ORFs between', orf_min, 'and', orf_max, 'in', fasta_filename,' = ', cnt