############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import itertools
# There exists pyfasta package -- http://pypi.python.org/pypi/pyfasta/
# Use it !

def get_lengths_from_fastafile(filename):
    """
        Gets filename of FASTA-file
        Returns list of lengths of sequences in FASTA-file
    """
    lengths = []
    l = 0
    for line in open(filename):
        if line[0] == '>':
            if l: # not the first sequence in FASTA
                lengths.append(l)
                l = 0
        else:
            l += len(line.strip())
    lengths.append(l)
    return lengths


def split_fasta(filename, outputdir):
    """
        Gets filename of FASTA-file and directory to output
        Creates separate FASTA-files for each sequence in FASTA-file
        Returns nothing
        Oops, similar to: pyfasta split --header "%(seqid).fasta" original.fasta
    """
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    outFile = None
    for line in open(filename):
        if line[0] == '>':
            if outFile:
                outFile.close()
            outFile = open(os.path.join(outputdir, line[1:].strip() + '.fa'), 'w')
        if outFile:
            outFile.write(line)
    if outFile: # if filename is empty
        outFile.close()


def read_fasta(filename):
    """
        Returns list of FASTA entries (in tuples: name, seq)
    """
    res_name = []
    res_seq = []
    first = True
    seq = ''
    fastafile = file
    file_ext = os.path.splitext(filename)[1]
    if file_ext == ".gz":
        import gzip
        fastafile = gzip.open(filename)
    else:
        fastafile = open(filename)

    for line in fastafile:
        if line[0] == '>':
            res_name.append(line.strip())
            if not first:
                res_seq.append(seq)
            else:
                first = False
            seq = ''
        else:
            seq += line.strip()
    res_seq.append(seq)
    return zip(res_name, res_seq)

def write_fasta(fasta):
    for name, seq in fasta:
        print name
        for i in xrange(0,len(seq),60):
            print seq[i:i+60]

def write_fasta_to_file(filename, fasta):
    outfile = open(filename, 'a')
    for name, seq in fasta:
        outfile.write(name + '\n')
        for i in xrange(0,len(seq),60):
            outfile.write(seq[i:i+60] + '\n')
    outfile.close()

def comp(letter):
    return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[letter.upper()]


def rev_comp(seq):
    return ''.join(itertools.imap(comp, seq[::-1]))

def remove_nonACGT(seq):
    seq2 = []
    for c in seq:
        if c  in 'ACGT':
	    seq2.append(c)
    return string.join(seq2, '')	
