import os
import gzip

def get_lengths_from_fastafile(filename):
    """
        Takes filename of FASTA-file
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
        Takes filename of FASTA-file and directory to output
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
    first = True
    seq = ''
    name = ''
    file_ext = os.path.splitext(filename)[1]
    fastafile = gzip.open(filename) if file_ext == ".gz" else open(filename)

    for line in fastafile:
        if line[0] == '>':
            if not first:
                yield name, seq
            first = False
            name = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
    if name or seq:
        yield name, seq

def write_fasta(fasta):
    for name, seq in fasta:
        print '>%s' % name
        for i in xrange(0,len(seq),60):
            print seq[i:i+60]

def write_fasta_to_file(filename, fasta):
    outfile = open(filename, 'w')
    for name, seq in fasta:
        outfile.write('>%s\n' % name)
        for i in xrange(0,len(seq),60):
            outfile.write(seq[i:i+60] + '\n')
    outfile.close()

